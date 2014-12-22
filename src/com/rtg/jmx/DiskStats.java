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

import static com.rtg.jmx.MonUtils.NF0;
import static com.rtg.jmx.MonUtils.NF2;

import java.io.IOException;

/**
 * Output disk utilization stats.
 */
public class DiskStats extends ExternalCommand implements MonStats {

  private static final String IO_STAT = "/usr/bin/iostat";
  private static final String NA = "n/a";
  private static final String HEADERBOT = " r/s  w/s  rMB/s  wMB/s  %util";

  private final String mDisk;

  private final boolean mEnabled;

  DiskStats(String disk) {
    super(IO_STAT, "-xk", "1", "2", disk);
    mDisk = disk;
    boolean enabled = false;
    try {
      final String[] result = runCommand(IO_STAT, "-xk", "1", "1", disk);
      if (result != null) {
        for (String s : result) {
          if (s.contains(mDisk)) {
            enabled = true;
          }
        }
      }
    } catch (IOException e) {
      // Ignore
    }
    mEnabled = enabled;
  }

  @Override
  public void addHeader(Appendable out) {  }

  @Override
  public void addColumnLabelsTop(Appendable out) throws IOException {
    if (!mEnabled) {
      return;
    }
    MonUtils.padRight(out, "---Disk-" + mDisk, HEADERBOT.length(), '-');
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
    String rs = NA;
    String ws = NA;
    String rms = NA;
    String wms = NA;
    String util = NA;
    if (lines != null) {
      for (int i = lines.length - 1; i >= 0; --i) {
        //System.err.println("DISK line " + lines[i]);
        if (lines[i].contains(mDisk)) {
          final String[] res = lines[i].split("[  ]+");
          rs = NF0.format(Double.valueOf(res[3]));
          ws = NF0.format(Double.valueOf(res[4]));
          rms = NF2.format(Double.valueOf(res[5]) / 1024);
          wms = NF2.format(Double.valueOf(res[6]) / 1024);
          util = NF2.format(Double.valueOf(res[11]));
          break;
        }
      }
    }
    final int width = 6;
    MonUtils.pad(out, rs, 4);
    out.append(" ");
    MonUtils.pad(out, ws, 4);
    out.append(" ");
    MonUtils.pad(out, rms, width);
    out.append(" ");
    MonUtils.pad(out, wms, width);
    out.append(" ");
    MonUtils.pad(out, util, width);
  }
}
