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

import java.io.Closeable;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Output stats from management beans.
 */
public class RecordStats implements Runnable {

  private static final String LS = System.lineSeparator();

  private ArrayList<MonStats> mStats = new ArrayList<>();

  private Appendable mOut;
  private int mDelay;
  private static final int REPEAT_HEADER = 30;
  private boolean mRun = true;


  RecordStats(Appendable out, int delay) {
    if (delay < 100) {
      throw new IllegalArgumentException();
    }
    mOut = out;
    mDelay = delay;
  }

  void addStats(MonStats... stats) {
    Collections.addAll(mStats, stats);
  }

  void terminate() {
    mRun = false;
  }

  void addHeader() throws IOException {
    mOut.append("# Monitoring started.").append(LS);
    for (MonStats s : mStats) {
      s.addHeader(mOut);
    }
  }

  void addColumnLabels() throws IOException {
    mOut.append("#");
    for (MonStats s : mStats) {
      mOut.append(" ");
      s.addColumnLabelsTop(mOut);
    }
    mOut.append(LS);
    mOut.append("#");
    for (MonStats s : mStats) {
      mOut.append(" ");
      s.addColumnLabelsBottom(mOut);
    }
    mOut.append(LS);
    mOut.append("#=============================================================================================================================").append(LS);
  }

  void addColumnData() throws IOException {
    mOut.append(" ");
    for (MonStats s : mStats) {
      mOut.append(" ");
      s.addColumnData(mOut);
    }
    mOut.append(LS);
  }

  @Override
  public void run() {
    try {
      addHeader();
      int count = 0;
      while (mRun) {
        final long time = System.currentTimeMillis();
        if (count == 0) {
          addColumnLabels();
        }
        addColumnData();

        count++;
        if (count == REPEAT_HEADER) {
          count = 0;
        }

        try {
          if (mRun) {
            final long wait = Math.max(100, mDelay - (System.currentTimeMillis() - time));
            Thread.sleep(wait);
          }
        } catch (InterruptedException e) {
          mRun = false;
        }
      }
      mOut.append("# Monitoring finished.").append(LS);
    } catch (IOException e) {
      System.err.println("Monitoring disabled: " + e.getMessage());
    } finally {
      if (mOut instanceof Closeable) {
        try {
          ((Closeable) mOut).close();
        } catch (IOException e) {
          // Do nothing
        }
      }
    }
  }

}
