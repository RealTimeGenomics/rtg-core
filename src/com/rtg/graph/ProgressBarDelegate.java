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

import javax.swing.JProgressBar;
import javax.swing.SwingUtilities;

import com.reeltwo.jumble.annotations.JumbleIgnore;

/**
 * Implementation that sends progress events to a progress bar
 */
@JumbleIgnore
public class ProgressBarDelegate implements ProgressDelegate {
  private final JProgressBar mProgressBar;

  private int mTotalLines;
  private int mTotalFiles;

  /**
   * @param prog progress bar to display progress on
   */
  public ProgressBarDelegate(JProgressBar prog) {
    this.mProgressBar = prog;
  }

  @Override
  public void setProgress(final int progress) {
    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        if (progress < 0) {
          mProgressBar.setIndeterminate(false);
        } else if (!mProgressBar.isStringPainted()) {
          mProgressBar.setStringPainted(true);
          mProgressBar.setIndeterminate(true);
          mProgressBar.setVisible(true);
        }
        mProgressBar.setString("" + progress);
      }
    });
  }

  @Override
  public void addFile(int numberLines) {
    mTotalFiles++;
    mTotalLines += numberLines;
  }

  @Override
  public void done() {
    setProgress(-1);
    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        mProgressBar.setString("Loaded " + mTotalLines + " points from " + mTotalFiles + " files.");
      }
    });
  }
}
