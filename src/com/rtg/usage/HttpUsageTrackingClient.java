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

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLEncoder;
import java.util.HashMap;
import java.util.Map;
import java.util.UUID;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.Environment;
import com.rtg.util.License;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.IOUtils;
import com.rtg.util.io.MemoryPrintStream;

/**
 * Client for talking to <code>http</code> based usage tracking server
 */
@TestClass(value = {"com.rtg.usage.UsageServerTest", "com.rtg.usage.HttpUsageTrackingClientTest"})
public class HttpUsageTrackingClient implements UsageTrackingClient {

  private static final int DEFAULT_RETRY_WAIT = 1000;
  private static final int NUM_TRIES = 3;
  private final String mHost;
  private final UsageConfiguration mUsageConfiguration;
  private final boolean mRequireUsage;

  /**
   * @param host <code>url</code> for posting usage to
   * @param conf the usage tracking configuration (for user and host name logging options)
   * @param requireUsage true if a failed message should be treated as an error
   */
  public HttpUsageTrackingClient(String host, UsageConfiguration conf, boolean requireUsage) {
    mHost = host;
    mUsageConfiguration = conf;
    mRequireUsage = requireUsage;
  }

  @Override
  public void recordEnd(long metric, String module, UUID runId, boolean success) {
    final HashMap<String, String> values = new HashMap<>();
    values.put(UsageServer.RUN_ID, runId.toString());
    values.put(UsageServer.MODULE, module);
    values.put(UsageServer.TYPE, success ? "Success" : "Fail");
    values.put(UsageServer.METRIC, Long.toString(metric));
    setCommon(values);
    sendMessage(values);
  }

  @Override
  public void recordBeginning(String module, UUID runId) {
    final HashMap<String, String> values = new HashMap<>();
    values.put(UsageServer.RUN_ID, runId.toString());
    values.put(UsageServer.TYPE, "Start");
    values.put(UsageServer.MODULE, module);
    //values.put(UsageServer.METRIC, "N/A");
    setCommon(values);
    sendMessage(values);
  }

  private void setCommon(Map<String, String> map) {
    if (mUsageConfiguration.logUsername()) {
      map.put(UsageServer.USERNAME, UsageMessage.trimField(System.getProperty("user.name"), UsageMessage.USERNAME_TRIM_LENGTH));
    }
    if (mUsageConfiguration.logHostname()) {
      map.put(UsageServer.HOSTNAME, UsageMessage.trimField(Environment.getHostName(), UsageMessage.HOSTNAME_TRIM_LENGTH));
    }
    if (mUsageConfiguration.logCommandLine()) {
      map.put(UsageServer.COMMANDLINE, UsageMessage.trimField(CommandLine.getCommandLine(), UsageMessage.COMMANDLINE_TRIM_LENGTH));
    }
    map.put(UsageServer.SERIAL, License.getSerialNumber());
    map.put(UsageServer.VERSION, Environment.getVersion());
  }

  HttpURLConnection openConnection() throws IOException {
    return (HttpURLConnection) new URL(mHost).openConnection();
  }

  int getRetryTimeInMillis() {
    return DEFAULT_RETRY_WAIT;
  }

  private void sendMessage(Map<String, String> values) {
    boolean success = false;
    String failureMessage = "";
    for (int t = 0; t < NUM_TRIES && !success; t++) {
      if (t > 0) {
        try {
          Thread.sleep(getRetryTimeInMillis());
        } catch (InterruptedException e) {
          if (mRequireUsage) {
            throw new NoTalkbackSlimException("Failed to send usage information. Aborting.");
          } else {
            return;
          }
        }
        //Diagnostic.warning("Retrying... (Attempt " + (t + 1) + " of " + NUM_TRIES + ")");
      }
      try {
        final HttpURLConnection http = openConnection();
        http.setConnectTimeout(30000);
        http.setReadTimeout(30000);
        http.setRequestMethod("POST");
        http.setDoOutput(true);
        http.setRequestProperty("Content-Type", "application/x-www-form-urlencoded");
        boolean first = true;
        final StringBuilder body = new StringBuilder();
        for (Map.Entry<String, String> entry : values.entrySet()) {
          if (entry.getValue() == null) {
            continue;
          }
          if (!first) {
            body.append("&");
          }
          first = false;
          body.append(URLEncoder.encode(entry.getKey(), "UTF-8"));
          body.append("=");
          body.append(URLEncoder.encode(entry.getValue(), "UTF-8"));
        }
        final int code;
        final String response;
        try (OutputStream os = http.getOutputStream()) {
          os.write(body.toString().getBytes());
        }
        code = http.getResponseCode();
        try (InputStream responseIn = http.getInputStream()) {
          response = IOUtils.readAll(responseIn);
        }
        if (code != 200 || !UsageServer.RTG_USAGE_ACCEPT.equals(response)) {
          failureMessage = "Failed to send usage information to " + mHost + " (" + code + " " + http.getResponseMessage() + ")";
          if (License.isDeveloper() || "Regression User".equals(License.getPerson())) {
            // For debugging failures during development / regression testing
            Diagnostic.warning("DEV:" + failureMessage);
          }
        } else {
          success = true; //all good!
        }
      } catch (IOException ioe) {
        failureMessage = "Failed to send usage information to " + mHost + " (" + ioe.getMessage() + ")";
        if (License.isDeveloper() || "Regression User".equals(License.getPerson())) {
          // For debugging failures during development / regression testing
          try (MemoryPrintStream mps = new MemoryPrintStream()) {
            mps.printStream().println(failureMessage);
            ioe.printStackTrace(mps.printStream());
            Diagnostic.warning("DEV:" + mps.toString());
          }
        }
      }
    }
    if (!success) {
      Diagnostic.warning(failureMessage);
      Diagnostic.warning("(Gave up after " + NUM_TRIES + " attempts)");
      if (mRequireUsage) {
        throw new NoTalkbackSlimException("Failed to send usage information. Aborting.");
      }
    }
  }
}
