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

import static com.rtg.jmx.MonUtils.NF2;

import java.io.IOException;

import com.rtg.util.Constants;

/**
 * Output network utilization stats.
 */
public class NetworkStats extends ExternalCommand implements MonStats {

  private static final String IFCONFIG = "/sbin/ifconfig";
  private static final String HEADERBOT = "txMB/s rxMB/s";
  private static final String NA = "n/a";
  private static final double MB = Constants.MB; // Convert to double


  private final boolean mEnabled;
  private final String mInterface;

  private long mTx = Long.MAX_VALUE;
  private long mRx = Long.MAX_VALUE;
  private long mTime = Long.MAX_VALUE;

  NetworkStats(String interfaceName) {
    super(IFCONFIG, interfaceName);
    boolean enabled = false;
    try {
      if (runCommand() != null) {
        enabled = true;
      }
    } catch (IOException e) {
      // Ignore
    }
    mEnabled = enabled;
    mInterface = interfaceName;
  }

  @Override
  public void addHeader(Appendable out) {  }

  @Override
  public void addColumnLabelsTop(Appendable out) throws IOException {
    if (!mEnabled) {
      return;
    }
    MonUtils.padRight(out, "--Net-" + mInterface, HEADERBOT.length(), '-');
  }

  @Override
  public void addColumnLabelsBottom(Appendable out) throws IOException {
    if (!mEnabled) {
      return;
    }
    out.append(HEADERBOT);
  }

  @Override
  public void addColumnData(Appendable out) throws IOException {
    if (!mEnabled) {
      return;
    }
    final String[] lines = runCommand();
    String txs = NA;
    String rxs = NA;
    if (lines != null) {
      for (int i = lines.length - 1; i >= 0; --i) {
        //System.err.println("network line : " + lines[i]);
        if (lines[i].contains("RX bytes")) {
          final String[] res = lines[i].trim().split("[  ]+");
          final long rx = Long.parseLong(res[1].split(":")[1]);
          final long tx = Long.parseLong(res[5].split(":")[1]);
          final long time = System.currentTimeMillis();
          final long dt = tx - mTx;
          final long dr = rx - mRx;
          final double dtime = (time - mTime) / 1000.0;
          mTx = tx;
          mRx = rx;
          mTime = time;
          if (dt >= 0 && dr >= 0) {
            txs = NF2.format(dt / dtime / MB);
            rxs = NF2.format(dr / dtime / MB);
          }
          break;
        }
      }
    }
    MonUtils.pad(out, txs, 6);
    out.append(" ");
    MonUtils.pad(out, rxs, 6);
  }

}
