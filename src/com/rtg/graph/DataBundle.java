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

import com.reeltwo.plot.Point2D;
import com.reeltwo.plot.PointPlot2D;
import com.reeltwo.plot.TextPlot2D;
import com.reeltwo.plot.TextPoint2D;

/**
 * Holds counts.
 */
final class DataBundle {
  private static final int TOTAL_LABELS = 10;

  private final Point2D[] mPoints;
  private final float[] mScores;

  private String mTitle;
  private final int mTotalVariants;
  private boolean mShow;

  private Point2D[] mRangedPoints = null;
  private float[] mRangedScores = null;
  private TextPoint2D[] mRangedPosPoints = null;

  DataBundle(String title, Point2D[] points, float[] scores, int totalVariants) {
    mPoints = points;
    mScores = scores;
    mTitle = title;
    mTotalVariants = totalVariants;
    mShow = true;

    resetRange();
  }

  boolean show() {
    return mShow;
  }

  void show(boolean flag) {
    mShow = flag;
  }

  void setTitle(String title) {
    mTitle = title;
  }

  void setScoreRange(float min, float max) {
    final int smin = (int) (min * mScores.length);
    final int smax = (int) (max * mScores.length);

    final ArrayList<Float> scores = new ArrayList<>();
    final ArrayList<Point2D> points = new ArrayList<>();

    for (int i = smin; i < smax; i++) {
      scores.add(mScores[i]);
      points.add(mPoints[i]);
    }

    mRangedPoints = points.toArray(new Point2D[points.size()]);
    mRangedScores = new float[scores.size()];
    for (int i = 0; i < scores.size(); i++) {
      mRangedScores[i] = scores.get(i);
    }

    updateLabels();
  }

  void resetRange() {
    mRangedPoints = mPoints;
    mRangedScores = mScores;
    updateLabels();
  }

  private void updateLabels() {
    final ArrayList<Integer> counts = new ArrayList<>();
    float px = 0;
    float py = 0;
    if (mRangedPoints.length != 0) {
      px = mRangedPoints[0].getX();
      py = mRangedPoints[0].getY();
    }
    int countTotal = 0;
    for (Point2D p : mRangedPoints) {
      final int c = (int) (p.getX() - px + p.getY() - py);
      counts.add(c);
      countTotal += c;
      px = p.getX();
      py = p.getY();
    }

    final ArrayList<TextPoint2D> posPoints = new ArrayList<>();
    // set up score labels - make TOTAL_LABELS per line
    if (countTotal != 0) {
      final int step = countTotal / TOTAL_LABELS;
      Point2D p = mRangedPoints[0];
      posPoints.add(new TextPoint2D(p.getX(), p.getY(), String.format("%.3g", mRangedScores[0])));
      if (step != 0) {
        int c = 0;
        for (int i = 0; i < counts.size(); i++) {
          for (int j = 0; j < counts.get(i); j++) {
            c++;
            if (c != 0 && c % step == 0) {
              p = mRangedPoints[i];
              posPoints.add(new TextPoint2D(p.getX(), p.getY(), String.format("%.3g", mRangedScores[i])));
            }
          }
        }
      }
      final int end = mRangedScores.length - 1;
      p = mRangedPoints[end];
      posPoints.add(new TextPoint2D(p.getX(), p.getY(), String.format("%.3g", mRangedScores[end])));
    }
    mRangedPosPoints = posPoints.toArray(new TextPoint2D[posPoints.size()]);
  }

  TextPoint2D getMaxRangedPoint() {
    return mRangedPosPoints.length == 0 ? null : mRangedPosPoints[mRangedPosPoints.length - 1];
  }

  PointPlot2D getPlot(int lineWidth, int colour) {
    final PointPlot2D lplot = new PointPlot2D();
    lplot.setData(mRangedPoints);
    lplot.setPoints(false);
    lplot.setLines(true);
    lplot.setLineWidth(lineWidth);
    lplot.setColor(colour);
    lplot.setTitle(mTitle);
    return lplot;
  }

  TextPlot2D getScoreLabels() {
    final TextPlot2D tplot = new TextPlot2D();
    tplot.setData(mRangedPosPoints);
    return tplot;
  }

  PointPlot2D getScorePoints(int lineWidth, int colour) {
    final PointPlot2D pplot = new PointPlot2D();
    pplot.setData(mRangedPosPoints);
    pplot.setPoints(true);
    pplot.setColor(colour);
    pplot.setLineWidth(lineWidth);
    return pplot;
  }

  int getTotalVariants() {
    return mTotalVariants;
  }

  public String getTitle() {
    return mTitle;
  }
}
