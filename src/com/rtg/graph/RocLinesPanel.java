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

import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.SwingUtilities;

import com.reeltwo.jumble.annotations.JumbleIgnore;

/**
 * A panel that holds a bunch of ROC line UI panels, handles changing the order of lines in the list.
 */
@JumbleIgnore
class RocLinesPanel extends Box {

  private final RocPlot mRocPlot;
  private final ArrayList<String> mPlotOrder = new ArrayList<>();

  RocLinesPanel(RocPlot rocPlot) {
    super(BoxLayout.Y_AXIS);
    mRocPlot = rocPlot;
  }

  ArrayList<String> plotOrder() {
    return mPlotOrder;
  }


  public void addLine(final RocLinePanel cp) {
    cp.addActionListener(new MoveActionListener(cp));
    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        mPlotOrder.add(cp.getText());
        add(cp);
      }
    });
  }

  @Override
  public Component add(Component comp) {
    final Component add = super.add(comp);
    int sizey = 0;
    for (Component c : getComponents()) {
      sizey += c.getPreferredSize().height;
    }
    setPreferredSize(new Dimension(getPreferredSize().width, sizey));
    return add;
  }

  @JumbleIgnore
  private class MoveActionListener implements ActionListener {
    private final RocLinePanel mPanel;

    public MoveActionListener(RocLinePanel panel) {
      mPanel = panel;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
      final Component[] components = getComponents();
      for (int i = 0; i < components.length; i++) {
        final RocLinePanel cp = (RocLinePanel) components[i];

        if (cp == mPanel) {
          if ("up".equals(e.getActionCommand())) {
            //System.err.println("Move up ");
            if (i > 0) {
              mPlotOrder.remove(i);
              mPlotOrder.add(i - 1, mPanel.getText());
              remove(mPanel);
              add(mPanel, i - 1);
            }
          } else {
            //System.err.println("Move down ");
            if (i < components.length - 1) {
              mPlotOrder.remove(i);
              mPlotOrder.add(i + 1, mPanel.getText());
              remove(mPanel);
              add(mPanel, i + 1);
            }
          }
          validate();
          mRocPlot.showCurrentGraph();
          break;
        }
      }
    }
  }

}
