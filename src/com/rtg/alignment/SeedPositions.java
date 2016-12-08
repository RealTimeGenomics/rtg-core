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
package com.rtg.alignment;

import com.rtg.mode.DnaUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * A helper class for storing seed positions and a few helper functions.
 *
 */
class SeedPositions {

  /** See ActionsHelper constants. */
  int mType;

  /** Start of read region, zero-based inclusive. */
  int mX1;

  /** End of read region, zero-based exclusive. */
  int mX2;

  /** Start of template region, zero-based inclusive. */
  int mY1;

  /** End of template region, zero-based exclusive. */
  int mY2;
  @Override
  public String toString() {
    return "[" + mX1 + " " + mX2 + ") [" + mY1 + " " + mY2 + ") ";
  }

  /**
   * @return the width of this region on the X axis.
   */
  protected int xWidth() {
    return mX2 - mX1;
  }

  /**
   * @return the width of this region on the Y axis.
   */
  protected int yWidth() {
    return mY2 - mY1;
  }

  /**
   * Dumps the seeds as a PPM file.
   *
   * @param seeds the list of seeds
   * @param f File
   * @param w width
   * @param h height
   * @throws IOException an IO exception if there is a problem
   */
/*  @JumbleIgnore
  public static void dumpSeedsAsPPM(final SeedPositions[] seeds, final File f, final int w, final int h) {
    final int[][][] ppm = makeArray(w + 1, h + 1, 3);

    for (int i = 0; i < seeds.length; ++i) {
      final int count = seeds[i].mX2 - seeds[i].mX1;
      int x = seeds[i].mX1;
      int y = seeds[i].mY1;
      final int t = seeds[i].mType;
      if (t == 1) {
        int j = 0;
        while ((j < count) && (x < w) && (y < h)) {
          if ((j == 0) || (j == (count - 1))) {
            // System.err.println(x + " " + y);
            ppm[x][y][0] = 255;
            ppm[x][y][1] = 0;
            ppm[x][y][2] = 0;
          } else {
            ppm[x][y][0] = 255;
            ppm[x][y][1] = 255;
            ppm[x][y][2] = 255;
          }
          ++x;
          ++y;
          ++j;
        }
      }
    }

    writePPM(f, w, h, ppm);
  }

  @JumbleIgnore
  private static int[][][] makeArray(final int w, final int h, final int z) {
    final int[][][] res = new int[w][][];
    for (int i = 0; i < w; ++i) {
      res[i] = new int[h][];
      for (int j = 0; j < h; ++j) {
        res[i][j] = new int[z];
      }
    }
    return res;
  }

  @JumbleIgnore
  private static void writePPM(final File f, final int w, final int h, final int[][][] ppm) {
    final PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f)));
    pw.print("P6 " + w + " " + h + " " + 255);
    pw.print("\r");
    for (int i = 0; i < w; ++i) {
      for (int j = 0; j < h; ++j) {
        final char r = (char) ppm[i][j][0];
        final char g = (char) ppm[i][j][1];
        final char b = (char) ppm[i][j][2];
        pw.print(r);
        pw.print(g);
        pw.print(b);
      }
    }
    pw.close();
  }
  */

  /**
   * Dumps the seeds to standard error
   *
   * @param seeds the seeds
   * @param numSeeds the number of seeds
   * @param read reads
   * @param template and the template
   */
  protected static void dumpSeeds(SeedPositions[] seeds, int numSeeds, byte[] read, byte[] template) {
    for (int i = 0; i < numSeeds; ++i) {
      System.err.println(seeds[i] + " " + DnaUtils.bytesToSequenceIncCG(read, seeds[i].mX1, seeds[i].mX2 - seeds[i].mX1)
          + " " + DnaUtils.bytesToSequenceIncCG(template, seeds[i].mY1, seeds[i].mY2 - seeds[i].mY1));
    }
  }

  /**
   * A pretty version of the seeds
   *
   * @param seeds the list of seeds
   * @return pretty string
   */
  protected static String dumpSeedsAsString(final SeedPositions[] seeds) {
    final StringBuilder sb = new StringBuilder();
    for (SeedPositions seed : seeds) {
      sb.append(seed.toString());
      sb.append(StringUtils.LS);
    }
    return sb.toString();
  }

  static void logBadSeed(String message, SeedPositions[] seeds, int i, byte[] read, int readlen, byte[] template, int zeroBasedStart) {
    Diagnostic.developerLog(message + StringUtils.LS + DnaUtils.bytesToSequenceIncCG(read) + StringUtils.LS + DnaUtils.bytesToSequenceIncCG(template, zeroBasedStart, readlen) + StringUtils.LS + DnaUtils.bytesToSequenceIncCG(read, seeds[i].mX1, seeds[i].xWidth()) + StringUtils.LS + DnaUtils.bytesToSequenceIncCG(template, zeroBasedStart + seeds[i].mY1, seeds[i].xWidth()) + StringUtils.LS);
  }

  static boolean seedIntegrity(final SeedPositions[] seeds, final int numSeeds, byte[] read, final int readlen, final byte[] template, final int zeroBasedStart) {
    // check ranges are increasing
    // System.err.println("seedintegrity");
    for (int i = 0; i < numSeeds; ++i) {
      if (seeds[i].mX1 > seeds[i].mX2) {
   //     SeedPositions.dumpSeeds(seeds, numSeeds, read, template);
        logBadSeed("SeededAligner failed integrity check: x -ve range", seeds, i, read, readlen, template, zeroBasedStart);
        return false;
      }
      if (seeds[i].mY2 - seeds[i].mY1 != seeds[i].mX2 - seeds[i].mX1) {
   //     SeedPositions.dumpSeeds(seeds, numSeeds, read, template);
        logBadSeed("SeededAligner failed integrity check: x and y ranges are different lengths", seeds, i, read, readlen, template, zeroBasedStart);
        return false;
      }
    }

    // check empty seed
    for (int i = 0; i < numSeeds; ++i) {
      if ((seeds[i].mX1 == seeds[i].mX2) && (seeds[i].mY1 == seeds[i].mY2)) {
        logBadSeed("SeededAligner failed integrity check: empty seed: " + seeds[i], seeds, i, read, readlen, template, zeroBasedStart);
        return false;
      }
    }

    // increase
    for (int i = 0; i < numSeeds - 1; ++i) {
      if (seeds[i + 1].mX1 <= seeds[i].mX1) {
        logBadSeed("SeededAligner failed integrity check: non monotonic x: " + seeds[i] + " " + seeds[i + 1] + " " + readlen + " " + zeroBasedStart, seeds, i, read, readlen, template, zeroBasedStart);
        return false;
      }
    }

    // also check off end of array
    for (int i = 0; i < numSeeds; ++i) {
      if ((seeds[i].mX2 > readlen) || (seeds[i].mY2 > template.length) || seeds[i].mX1 < 0 || seeds[i].mY1 < 0) {
        logBadSeed("SeededAligner failed integrity check: BUG: off end of array", seeds, i, read, readlen, template, zeroBasedStart);
        return false;
      }
    }

    // check the seed is really equal
    for (int i = 0; i < numSeeds; ++i) {
      int start = seeds[i].mY1;
      for (int j = seeds[i].mX1; j < seeds[i].mX2; ++j) {
        if (read[j] != template[start]) {
          logBadSeed("SeededAligner failed integrity check: bad seed at " + j + " " + start + "   values = " + read[j] + " " + template[start] + " " + seeds[i] + " " + StringUtils.LS, seeds, i, read, readlen, template, zeroBasedStart);
          return false;
        }
        ++start;
      }
    }
    return true;
  }

  static boolean overlapIntegrity(SeedPositions[] seeds, int numSeeds, byte[] read, int readlen, byte[] template, int zeroBasedStart) {
    // don't overlap with the next seed
    for (int i = 0; i < numSeeds - 1; ++i) {
      if ((seeds[i].mX2 > seeds[i + 1].mX1) || (seeds[i].mY2 > seeds[i + 1].mY1)) {
//        SeedPositions.dumpSeeds(seeds, numSeeds, read, template);
        logBadSeed("SeededAligner failed integrity check: overlapping seeds", seeds, i, read, readlen, template, zeroBasedStart);
        return false;
      }
    }
    return true;
  }
}
