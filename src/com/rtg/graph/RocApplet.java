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

import java.awt.BorderLayout;
import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

import javax.swing.JApplet;

import com.reeltwo.jumble.annotations.JumbleIgnore;
import com.rtg.util.gzip.WorkingGzipInputStream;
import com.rtg.util.io.FileUtils;

/**
 * An applet version of the ROC plot.
 */
@JumbleIgnore
public class RocApplet extends JApplet {

  /* Embed in HTML with something like:

<applet width="100%" height="100%" code="com.rtg.graph.RocApplet" archive="ROC.jar">
<param name="data1" value="snp_pe_snp1_cal.eval/weighted_roc.tsv"> <param name="name1" value="curve number 1">
<param name="data2" value="snp_pe_snp2_cal.eval/weighted_roc.tsv"> <param name="name2" value="funky curve">
<param name="data3" value="snp_pe_snp2_nocal.eval/weighted_roc.tsv"> <param name="name3" value="another curve">
<param name="data4" value="snp_pe_snp1_nocal.eval/weighted_roc.tsv">
</applet>

   */

  static final String PARAM_SIDEPANE = "sidepane";
  static final String PARAM_NAME = "name";
  static final String PARAM_DATA = "data";
  static final String PARAM_TITLE = "title";
  static final String PARAM_SCORES = "scores";
  static final String PARAM_LINEWIDTH = "linewidth";

  RocPlot mRocPlot = null;

  @Override
  public void start() {
    try {
      javax.swing.SwingUtilities.invokeAndWait(new Runnable() {
        @Override
        public void run() {
          mRocPlot.showCurrentGraph();
        }
      });
    } catch (final Exception e) {
      System.err.println("createGUI didn't successfully complete: " + e.getMessage());
    }
  }

  @Override
  public void init() {
    try {
      mRocPlot = new RocPlot();
      javax.swing.SwingUtilities.invokeLater(new Runnable() {
        @Override
        public void run() {
          setLayout(new BorderLayout());
          add(mRocPlot.getMainPanel(), BorderLayout.CENTER);
          setGlassPane(mRocPlot.getZoomPlotPanel());
          getGlassPane().setVisible(true);
          mRocPlot.showOpenButton(false);

          final String title = getParameter(PARAM_TITLE);
          if (title != null) {
            System.err.println("Setting title to " + title);
            mRocPlot.setTitle(title);
          }

          final boolean scores = Boolean.parseBoolean(getParameter(PARAM_SCORES));
          System.err.println("Setting show scores to " + scores);
          mRocPlot.showScores(scores);

          final boolean sidepane = Boolean.parseBoolean(getParameter(PARAM_SIDEPANE));
          System.err.println("Setting sidepane visibility to " + sidepane);
          if (sidepane) {
            mRocPlot.setSplitPaneDividerLocation(1.0);
          }

          final String lineWidth = getParameter(PARAM_LINEWIDTH);
          if (lineWidth != null) {
            System.err.println("Setting line width to " + lineWidth);
            mRocPlot.setLineWidth(Integer.parseInt(lineWidth));
          }
}
      });

      final Thread t = new Thread() {
        @Override
        public void run() {
          // Get parameters and set them loading
          for (int streamNum = 1; getParameter(PARAM_DATA + streamNum) != null; streamNum++) {
            final String source = getParameter(PARAM_DATA + streamNum);
            String name = getParameter(PARAM_NAME + streamNum);
            if (name == null) {
              name = source;
            }
            try {
              System.err.println("Attempting to load data source: " + source);
              final URL sourceUrl = new URL(getCodeBase(), source);
              System.err.println("SourceUrl: " + sourceUrl);
              InputStream is = sourceUrl.openStream();
              if (FileUtils.isGzipFilename(source)) {
                is = new WorkingGzipInputStream(is);
              }
              ParseRocFile.loadStream(mRocPlot.getProgressBarDelegate(), new BufferedInputStream(is), name, true);
            } catch (final IOException e) {
              final String message = "Could not load data source: " + source + " " + e.getMessage();
              System.err.println(message);
              // TODO set status in UI
              mRocPlot.setStatus(message);
            }
          }
          // give gui a chance to catch up...
          try {
            Thread.sleep(1000);
          } catch (InterruptedException ie) {
          }
          mRocPlot.updateProgress();
          mRocPlot.showCurrentGraph();
        }
      };
      t.start();
    } catch (Exception e) {
      System.err.println("createGUI didn't successfully complete: " + e.getMessage());
    }
  }


}
