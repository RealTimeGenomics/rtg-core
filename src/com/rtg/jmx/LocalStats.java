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
package com.rtg.jmx;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;

/**
 * Set up a separate stats monitoring thread from within the same JVM
 * as SLIM.
 *
 */
public final class LocalStats {

  /** Name of property specifying monitor output destination */
  public static final String MON_DEST = "rtg.jmxmon";
  /** Destination specifier for within output directory */
  public static final String MON_DEST_OUTDIR = "outdir";

  /** Delay in seconds between updates */
  private static final int MON_DELAY = Integer.parseInt(System.getProperty("rtg.jmxmon.delay", "10"));

  /** Comma separated list of hard disks to monitor */
  private static final String MON_DISK = System.getProperty("rtg.jmxmon.disk", "sda");

  /** Comma separated list of network interfaces to monitor */
  private static final String MON_NET = System.getProperty("rtg.jmxmon.net", "eth0");


  private static RecordStats sStats = null;
  private static Thread sThread = null;


  private LocalStats() { }


  // Create a RecordStats object using values defined in system properties
  static RecordStats getLocalStats() {
    return getLocalStats(System.getProperty(MON_DEST), MON_DELAY, MON_DISK, MON_NET);
  }

  /**
   * Make a RecordStats pointing at suitable local object.
   *
   * @param dest destination, either <code>auto</code>, <code>err</code>, or a file name
   * @param delay delay in seconds.
   * @param disk comma separated list of disk names to monitor
   * @param net comma separated list of network interface names to monitor
   * @return a <code>RecordStats</code> value
   */
  static RecordStats getLocalStats(String dest, int delay, String disk, String net) {
    RecordStats rs = null;
    if (dest != null) {
      try {
        final PrintStream monout;
        if (dest.equals("auto")) {
          monout = new PrintStream(new FileOutputStream(File.createTempFile("jmxmon", ".txt", new File(System.getProperty("user.dir")))), true);
        } else if (dest.equals("err")) {
          monout = System.err;
        } else if (dest.equals("outdir")) {
          System.err.println("Output directory has not been set");
          throw new IOException();
        } else {
          monout = new PrintStream(new File(dest));
        }
        rs = new RecordStats(monout, delay * 1000);
        rs.addStats(new MBeanStats());
        for (String d : disk.split(",")) {
          rs.addStats(new DiskStats(d));
        }
        for (String n : net.split(",")) {
          rs.addStats(new NetworkStats(n));
        }
        rs.addStats(new ProgressStats());
      } catch (IOException e) {
        // Do nothing
      }
    }
    return rs;
  }

  /**
   * Starts monitoring in a separate thread. If <code>MON_DEST</code>
   * is not set or is set to <code>MON_DEST_OUTDIR</code>, or if
   * monitoring has already been initiated, this method does nothing.
   */
  public static synchronized void startRecording() {
    if (sThread != null) {
      return;
    }
    final String dest = System.getProperty(MON_DEST);
    if ((dest == null) || dest.equals(MON_DEST_OUTDIR)) {
      return;
    }
    sThread = new Thread(new Runnable() {
        @Override
        public void run() {
          sStats = getLocalStats();
          if (sStats != null) {
            System.err.println("Starting monitoring");
            sStats.run();
          }
        }
      });
    sThread.setDaemon(true);
    sThread.start();
  }

  /**
   * Asynchronously instructs the monitoring thread to stop.
   */
  public static synchronized void stopRecording() {
    if (sStats != null) {
      System.err.println("Shutting down monitoring");
      sStats.terminate();
      sThread.interrupt();
      Thread.yield(); // Give the other one more time to shut down.
    }
  }

  /**
   * Test out the local stats monitoring.
   *
   * @param args ignored
   * @exception IOException if an error occurs.
   * @exception InterruptedException if an error occurs.
   */
  public static void main(String[] args) throws IOException, InterruptedException {
    final RecordStats rs = getLocalStats("err", 5, "sda", "eth0");
    rs.run();
  }

}
