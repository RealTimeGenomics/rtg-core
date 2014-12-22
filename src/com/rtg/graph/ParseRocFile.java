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
package com.rtg.graph;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

import com.reeltwo.plot.Point2D;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 */
public final class ParseRocFile {

  private ParseRocFile() {
  }

  /**
   * loads ROC file into data bundle
   * @param progressBarDelegate called every 100 lines with progress, and at end with file stats
   * @param is input data. this stream is closed by this method.
   * @param shortName name for line
   * @param showProgress true if progress should be sent
   * @return the data bundle
   * @throws IOException if an IO error occurs
   */
  static DataBundle loadStream(ProgressDelegate progressBarDelegate, final BufferedInputStream is, final String shortName, boolean showProgress) throws IOException {
    int lines = 0;
    int totalVariants = -1;
    final ArrayList<Point2D> points = new ArrayList<>();
    final ArrayList<Float> scores = new ArrayList<>();

    String line = null;
    try (BufferedReader br = new BufferedReader(new InputStreamReader(is))) {
      float prevScore = Float.MIN_VALUE;
      float prevX = 0.0f;
      float prevY = 0.0f;
      float score = 0.0f;
      while ((line = br.readLine()) != null) {
        if (line.startsWith("#")) {
          if (line.contains("#total")) {
            final String[] parts = line.split("\\s");
            totalVariants = Integer.parseInt(parts[parts.length - 1]);
          }
          continue;
        }
        final String[] linearr = line.split("\t");
        if (linearr.length < 3) {
          throw new NoTalkbackSlimException("Malformed line: " + line + " in \"" + shortName + "\"");
        }
        final float x = Float.parseFloat(linearr[2]); // False positives
        final float y = Float.parseFloat(linearr[1]); // True positives
        score = Float.parseFloat(linearr[0]); // Score
        if (Float.compare(score, prevScore) != 0) {
          points.add(new Point2D(prevX, prevY));
          scores.add(score);
        }
        prevX = Math.max(prevX, x);
        prevY = Math.max(prevY, y);
        prevScore = score;

        lines++;
        if (showProgress && lines % 100 == 0) {
          progressBarDelegate.setProgress(lines);
        }
      }
      points.add(new Point2D(prevX, prevY));
      scores.add(score);
    } catch (final NumberFormatException e) {
      throw new NoTalkbackSlimException("Malformed line: " + line + " in \"" + shortName + "\"");
    }
    final float[] scoresArr = new float[scores.size()];
    for (int i = 0; i < scores.size(); i++) {
      scoresArr[i] = scores.get(i);
    }
    progressBarDelegate.addFile(lines);
    return new DataBundle(shortName, points.toArray(new Point2D[points.size()]), scoresArr, totalVariants);
  }
}
