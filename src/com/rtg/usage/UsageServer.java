/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.usage;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.net.InetSocketAddress;
import java.net.URLDecoder;
import java.nio.channels.FileLock;
import java.nio.charset.Charset;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.IOUtils;
import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;

/**
 * Server for <code>http</code> based usage tracking
 */
public class UsageServer {

  private final Object mSync = new Object();

  static final String RTG_USAGE_ACCEPT = "RTG Usage Accept";
  static final String RTG_USAGE_REJECT = "RTG Usage Reject";

  static final String SERIAL = "serial";
  static final String RUN_ID = "run_id";
  static final String VERSION = "version";
  static final String MODULE = "module";
  static final String METRIC = "metric";
  static final String TYPE = "type";
  static final String USERNAME = "username";
  static final String HOSTNAME = "hostname";
  static final String COMMANDLINE = "cmdline";
  private static final String UTF_8 = "UTF-8";

  private final File mUsageDir;
  private final int mPort;
  private final ThreadPoolExecutor mThreadPoolExecutor;

  private HttpServer mServer = null;
  private RandomAccessFile mRandomAccessFile;
  private FileLock mLock;
  private String mLastKey;
  private File mCurrentUsageFile;

  private int mStopTimer = 5;


  /**
   * Constructs and obtains lock on usage file for current month.
   *
   * @param port port to use for connections
   * @param usageDir directory to write usage file to
   * @param threads number of concurrent connections
   * @throws IOException if it happens
   */
  public UsageServer(int port, File usageDir, int threads) throws IOException {
    if (!Charset.isSupported(UTF_8)) {
      //in theory this should never ever happen. UTF-8 is a standard required charset by the java platform
      throw new IOException("UTF-8 character encoding not supported");
    }
    //flaky stuff
    System.setProperty("sun.net.httpserver.maxReqTime", "30");
    System.setProperty("sun.net.httpserver.maxRspTime", "30");
    mUsageDir = usageDir;
    mPort = port;
    mThreadPoolExecutor = new ThreadPoolExecutor(threads, threads, 0, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>(), new UsageServerThreadFactory());
  }

  /**
   * start the server
   * @throws IOException if it happens
   */
  public void start() throws IOException {
    switchToFile(FileUsageTrackingClient.ensureUsageFile(mUsageDir, getDate()));
    try {
      mServer = HttpServer.create(new InetSocketAddress(mPort), 50);
    } catch (IOException e) {
      throw new NoTalkbackSlimException("Could not create socket on port " + mPort + ". " + e.getMessage());
    }
    mServer.createContext("/", new UsageHandler());
    mServer.setExecutor(mThreadPoolExecutor);
//    mLastKey = FileUsageTrackingClient.getLastKey(mRandomAccessFile);
//    mRandomAccessFile.seek(mRandomAccessFile.length());
    mServer.start();
  }

  /**
   * @return obtain the
   * port in use by the server
   */
  public int getPort() {
    return mServer.getAddress().getPort();
  }

  /**
   * stop the server and release any file locks
   * @throws IOException if it happens
   */
  public void end() throws IOException {
    mServer.stop(mStopTimer);
    synchronized (mSync) {
      closeCurrentFile();
    }
    mThreadPoolExecutor.shutdown();
  }

  void setStopTimer(int stopTimer) {
    mStopTimer = stopTimer;
  }

  private void closeCurrentFile() throws IOException {
    if (mLock != null) {
      mLock.release();
      mRandomAccessFile.close();
      mLock = null;
      mRandomAccessFile = null;
      mLastKey = null;
    }
  }

  private void switchToFile(File fileToUse) throws IOException {
    closeCurrentFile();
    mCurrentUsageFile = fileToUse;
    mRandomAccessFile = new RandomAccessFile(mCurrentUsageFile, "rw");
    mLock = mRandomAccessFile.getChannel().tryLock();
    if (mLock == null) {
      throw new IOException("Could not lock usage file as it is already locked by another process");
    }
    mLastKey = FileUsageTrackingClient.getLastKey(mRandomAccessFile);
    mRandomAccessFile.seek(mRandomAccessFile.length());
  }

  private static class UsageServerThreadFactory implements ThreadFactory {
    private final AtomicInteger mThreadNo = new AtomicInteger(0);

    @Override
    public Thread newThread(Runnable r) {
      final Thread t = new Thread(r, "UsageServerThread-" + mThreadNo.getAndIncrement());
      if (t.isDaemon()) {
        t.setDaemon(false);
      }
      return t;
    }
  }

  private class UsageHandler implements HttpHandler {

    @Override
    public void handle(HttpExchange httpExchange) throws IOException {
      //post items exist in plain text in the request body
      final String postString;
      try (final InputStream requestBody = httpExchange.getRequestBody()) {
        postString = IOUtils.readAll(requestBody);
      }
      final Map<String, String> post;
      if ("POST".equals(httpExchange.getRequestMethod())) {
        post = parsePost(postString);
      } else {
        post = new HashMap<>();
      }
      final boolean messageOk = post.containsKey(SERIAL) && post.containsKey(RUN_ID) && post.containsKey(VERSION) && post.containsKey(MODULE) && post.containsKey(TYPE);
      final String reply = messageOk ? RTG_USAGE_ACCEPT : RTG_USAGE_REJECT;
      final Date d = getDate(); //before response so start/end calls an be guaranteed to be in ascending order
      httpExchange.sendResponseHeaders(200, reply.getBytes().length);
      try (OutputStream out = httpExchange.getResponseBody()) {
        out.write(reply.getBytes());
      }
      httpExchange.close();
      final UsageMessage message = UsageMessage.setMessage(post.get(SERIAL), post.get(RUN_ID), post.get(VERSION), post.get(MODULE), post.containsKey(METRIC) ? post.get(METRIC) : null, post.get(TYPE));
      if (post.containsKey(USERNAME)) {
        message.setUsername(post.get(USERNAME));
      }
      if (post.containsKey(HOSTNAME)) {
        message.setHostname(post.get(HOSTNAME));
      }
      if (post.containsKey(COMMANDLINE)) {
        message.setCommandLine(post.get(COMMANDLINE));
      }
      synchronized (mSync) {
        final File fileToUse = FileUsageTrackingClient.ensureUsageFile(mUsageDir, d);
        if (!fileToUse.equals(mCurrentUsageFile)) {
          switchToFile(fileToUse);
        }
        recordMessage(d, message);
      }
    }
  }

  /**
   * This gets called in the constructor. Please be careful if overriding (i.e. please don't override outside of tests)
   * @return the current date
   */
  Date getDate() {
    return new Date();
  }

  private void recordMessage(Date d, UsageMessage message) throws IOException {
    message.setDate(d);
    mRandomAccessFile.write(message.formatLine(mLastKey).getBytes());
    mLastKey = message.getChecksum();
  }

  private static Map<String, String> parsePost(String post) throws IOException {
    final HashMap<String, String> ret = new HashMap<>();
    if (post.length() > 0) {
      final String[] vals = StringUtils.split(post, '&');
      for (String val : vals) {
        final String[] keyval = StringUtils.split(val, '=');
        if (keyval.length == 2) {
          try {
            ret.put(URLDecoder.decode(keyval[0], UTF_8), URLDecoder.decode(keyval[1], UTF_8));
          } catch (IllegalArgumentException e) {
            Diagnostic.userLog("Failed to decode parameter: '" + val + "' (" + e.getMessage() + ")");
          }
        }
      }
    }
    return ret;
  }
}
