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
import java.io.RandomAccessFile;
import java.nio.channels.FileLock;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.UUID;

import com.rtg.util.Environment;
import com.rtg.util.License;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Implements a file based usage tracking client and server. Files written are based on the month and year
 */
public class FileUsageTrackingClient implements UsageTrackingClient {

  private static final int BUF_SIZE = 1024;
  private final File mUsageDir;
  private final UsageConfiguration mUsageConfiguration;
  private final boolean mRequireUsage;

  /**
   * @param usageDir directory in which to write usage tracking files
   * @param conf the usage tracking configuration (for user and host name logging options)
   * @param requireUsage true if usage tracking is mandatory, affects how error conditions are handled.
   */
  public FileUsageTrackingClient(File usageDir, UsageConfiguration conf, boolean requireUsage) {
    mUsageDir = usageDir;
    mUsageConfiguration = conf;
    mRequireUsage = requireUsage;
  }

  @Override
  public void recordBeginning(String module, UUID runId) {
    recordMessage(UsageMessage.startMessage(License.getSerialNumber(), runId.toString(), Environment.getVersion(), module));
  }

  @Override
  public void recordEnd(long metric, String module, UUID runId, boolean success) {
    recordMessage(UsageMessage.setMessage(License.getSerialNumber(), runId.toString(), Environment.getVersion(), module, Long.toString(metric), success ? "Success" : "Fail"));
  }

  private void recordMessage(UsageMessage message) {
    final Date d = new Date();
    message.setDate(d);
    if (mUsageConfiguration.logHostname()) {
      message.setHostname(Environment.getHostName());
    }
    if (mUsageConfiguration.logUsername()) {
      message.setUsername(System.getProperty("user.name"));
    }
    if (mUsageConfiguration.logCommandLine()) {
      message.setCommandLine(CommandLine.getCommandLine());
    }
    try {
      final File usageOut = ensureUsageFile(mUsageDir, d);
      try (final RandomAccessFile raf = new RandomAccessFile(usageOut, "rw");
           final FileLock lock = raf.getChannel().lock()) { //this blocks until lock is available, no timeout setting
        //get previous signage key thingy
        final String prevKey = getLastKey(raf);
        //write message
        raf.seek(raf.length());
        raf.write(message.formatLine(prevKey).getBytes());
        lock.release(); //this is because compiler complains about the lack of usage on this variable. Additionally it makes the intent clear about when this lock should be released
      }
    } catch (IOException e) {
      final String errorMessage = "Failed to record usage information (" + e.getMessage() + ")";
      if (mRequireUsage) {
        throw new NoTalkbackSlimException(errorMessage);
      } else {
        Diagnostic.warning(errorMessage);
      }
    }
  }

  static File ensureUsageFile(File usageDir, Date d) throws IOException {
    final File usageFile = getUsageFile(usageDir, d);
    if (!usageFile.exists()) {
      if (!usageDir.exists()) {
        if (!usageDir.mkdir()) {
          throw new IOException("Cannot create usage directory: " + usageDir.toString());
        } else {
          if (!usageDir.setReadable(true, false)
              || !usageDir.setWritable(true, false)) {
            throw new IOException("Failed to set permissions on usage directory:" + usageDir.toString());
          }
        }
      }
      // create file and set permissions so any user can write to it
      // file is created during random access and lock process
      try (final RandomAccessFile raf = new RandomAccessFile(usageFile, "rw");
          final FileLock lock = raf.getChannel().lock()) {

        if (!usageFile.setReadable(true, false)
            || !usageFile.setWritable(true, false)) {
          throw new IOException("Failed to set permissions on usage file:" + usageFile.toString());
        }
        lock.release();
      }
    }
    return usageFile;
  }

  static File getUsageFile(File usageDir, Date d) {
    final DateFormat fileDateFormat = new SimpleDateFormat("yyyy-MM");
    return new File(usageDir, fileDateFormat.format(d) + ".usage");
  }

  static String getLastKey(RandomAccessFile raf) throws IOException {
    return getLastKey(raf, BUF_SIZE);
  }
  static String getLastKey(RandomAccessFile raf, int bufSize) throws IOException {
    final StringBuilder sb = new StringBuilder();
    final byte[] buf = new byte[bufSize];
    long endPos = raf.length();
    while (endPos != 0) {
      final long seekPos = Math.max(endPos - buf.length, 0);
      raf.seek(seekPos);
      final int length = (int) (endPos - seekPos);
      raf.readFully(buf, 0, length);
      sb.insert(0, new String(buf, 0, length));
      final String currentString = sb.toString();
      final int index = currentString.lastIndexOf(UsageMessage.SIGNATURE);
      if (index == -1) {
        final int chopIndex = Math.max(currentString.indexOf('\r'), currentString.indexOf('\n'));
        sb.delete(chopIndex + 1, sb.length());
      } else {
        final int lastSlashN = currentString.indexOf('\n', index);
        final int lastSlashR = currentString.indexOf('\r', index);
        //trying to find everything from the start of the signature up to the start of the newline (which can be \r\n, \n or \r)
        //not found == -1, so if only one was -1 take the larger (i.e. found) one, if both were found take the smaller (first found) one. If both are -1 it doesn't matter.
        final int lastNewline = lastSlashR != -1 && lastSlashN != -1 ? Math.min(lastSlashR, lastSlashN) : Math.max(lastSlashR, lastSlashN);
        final int endSig = lastNewline != -1 ? lastNewline : currentString.length();
        return currentString.substring(index + UsageMessage.SIGNATURE.length(), endSig);
      }
      endPos = seekPos;
    }
    return null;
  }

}
