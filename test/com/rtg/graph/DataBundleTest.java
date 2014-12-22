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

import java.util.ArrayList;

import com.reeltwo.plot.Datum2D;
import com.reeltwo.plot.Graph2D;
import com.reeltwo.plot.Point2D;
import com.reeltwo.plot.PointPlot2D;
import com.reeltwo.plot.TextPoint2D;

import junit.framework.TestCase;

/**
 */
public class DataBundleTest extends TestCase {
  public void test() {
    final Point2D[] points = {new Point2D(4.0f, 5.0f), new Point2D(100.0f, 200.0f)};
    final float[] scores = {5.0f, 9.0f};
    final DataBundle db = new DataBundle("Monkey", points, scores, 300);
    assertEquals(300, db.getTotalVariants());
    assertEquals(100.0, db.getPlot(1, 1).getHi(Graph2D.X), 1e-9);
    assertEquals(200.0, db.getPlot(1, 1).getHi(Graph2D.Y), 1e-9);
    assertEquals(4.0, db.getPlot(1, 1).getLo(Graph2D.X), 1e-9);
    assertEquals(5.0, db.getPlot(1, 1).getLo(Graph2D.Y), 1e-9);
    assertEquals("Monkey", db.getTitle());
    assertEquals(2, db.getPlot(1, 1).getData().length);
  }

  public void testScoreLabels() {
    final ArrayList<Point2D> points = new ArrayList<>();
    final ArrayList<Float> scores = new ArrayList<>();
    for (int i = 0; i < 100; i++) {
      points.add(new Point2D(i, i));
      scores.add((float) (100 - i));
    }
    final float[] scoresArr = new float[scores.size()];
    for (int i = 0; i < scoresArr.length; i++) {
      scoresArr[i] = scores.get(i);
    }
    final DataBundle db = new DataBundle("Monkey", points.toArray(new Point2D[points.size()]), scoresArr, 300);
    db.setScoreRange(0.0f, 1.0f);
    final PointPlot2D scorePoints = db.getScorePoints(1, 1);
    final String[] exp = {"100", "90.0", "81.0", "71.0", "62.0", "52.0", "43.0", "33.0", "24.0", "14.0", "5.00", "1.00"};
    final Datum2D[] data = scorePoints.getData();
    for (int i = 0; i < data.length; i++) {
      final Datum2D d = data[i];
      assertTrue(d instanceof TextPoint2D);
      final TextPoint2D p = (TextPoint2D) d;
      assertEquals(exp[i], p.getText());
    }
  }
}
