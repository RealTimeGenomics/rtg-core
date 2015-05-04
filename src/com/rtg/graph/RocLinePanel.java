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

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JCheckBox;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.plaf.basic.BasicArrowButton;

import com.reeltwo.jumble.annotations.JumbleIgnore;
import com.reeltwo.plot.TextPoint2D;

/**
 * A panel that lets you enable or disable showing of a ROC line, rename the ROC line, and has buttons for
 * moving the line up or down in the list.
 */
@JumbleIgnore
class RocLinePanel extends JPanel {

  private final JCheckBox mCheckBox;
  private final BasicArrowButton mUpButton;
  private final BasicArrowButton mDownButton;
  private final DataBundle mDataBundle;
  private final String mPath;
  private final JTextField mTextField;
  private final RocPlot mRocPlot;
  private final JProgressBar mStatusBar;

  RocLinePanel(RocPlot rocPlot, final String path, final String name, DataBundle data, JProgressBar statusBar) {
    super(new GridBagLayout());
    mRocPlot = rocPlot;
    mPath = path;
    mDataBundle = data;
    mCheckBox = new JCheckBox();
    mCheckBox.setSelected(true);
    mTextField = new JTextField();
    mTextField.setText(name);
    mTextField.setEditable(false);
    mStatusBar = statusBar;
    mTextField.addMouseMotionListener(new MouseMotionAdapter() {
      @Override
      public void mouseMoved(MouseEvent e) {
        mStatusBar.setString(getTooltip(mPath, mTextField.getText()));
      }
    });

    mTextField.addFocusListener(new FocusAdapter() {

      @Override
      public void focusLost(FocusEvent e) {
        //mData.get(path).setTitle(mTextField.getText());
        mTextField.setText(mDataBundle.getTitle());
        mTextField.setEditable(false);
        mTextField.getCaret().setVisible(false);
        mRocPlot.showCurrentGraph();
      }
    });


    mTextField.addKeyListener(new KeyAdapter() {

      @Override
      public void keyPressed(KeyEvent e) {
        if (e.getKeyCode() == KeyEvent.VK_ESCAPE) {
          mTextField.setText(mDataBundle.getTitle());
          mTextField.setEditable(false);
          mTextField.getCaret().setVisible(false);
        }
        if (e.getKeyCode() == KeyEvent.VK_ENTER) {
          mDataBundle.setTitle(mTextField.getText());
          mTextField.setEditable(false);
          mTextField.getCaret().setVisible(false);
          mRocPlot.showCurrentGraph();
        }
      }

    });

    mTextField.addMouseListener(new MouseAdapter() {
      @Override
      public void mouseClicked(MouseEvent e) {
        if (e.getClickCount() == 2) {
          mTextField.getCaret().setVisible(true);
          mTextField.setEditable(true);
          mTextField.requestFocus();
        } else if (e.getClickCount() == 1 && !mTextField.isEditable()) {
          setSelected(!mCheckBox.isSelected());
        }
      }
    });

    final Box textPanel = Box.createHorizontalBox();
    textPanel.add(mCheckBox);
    textPanel.add(mTextField);
    final GridBagConstraints textConstraints = new GridBagConstraints();
    textConstraints.gridx = 0; textConstraints.gridy = 0;
    textConstraints.fill = GridBagConstraints.HORIZONTAL; textConstraints.weightx = 2;
    add(textPanel, textConstraints);

    final JPanel updownPanel = new JPanel();
    new BoxLayout(updownPanel, BoxLayout.X_AXIS);

    mUpButton = new BasicArrowButton(BasicArrowButton.NORTH);
    mUpButton.setActionCommand("up");
    mDownButton = new BasicArrowButton(BasicArrowButton.SOUTH);
    mDownButton.setActionCommand("down");

    updownPanel.add(mUpButton);
    updownPanel.add(mDownButton);
    updownPanel.setMinimumSize(updownPanel.getPreferredSize());
    final GridBagConstraints updownConstraints = new GridBagConstraints();
    updownConstraints.gridx = 1; updownConstraints.gridy = 0;
//    updownConstraints.weightx = 2;
    add(updownPanel, updownConstraints);
    mCheckBox.addItemListener(new CheckBoxListener(this));

    final JSlider rangeSlider = new JSlider();
//    rangeSlider.setPreferredSize(new Dimension(240, rangeSlider.getPreferredSize().height));
    rangeSlider.setMinimum(0);
    rangeSlider.setMaximum(1000);
    rangeSlider.setValue(1000);

    rangeSlider.addChangeListener(new ChangeListener() {
      @Override
      public void stateChanged(ChangeEvent e) {
        final JSlider slider = (JSlider) e.getSource();
        mDataBundle.setScoreRange(0.0f, slider.getValue() / 1000.0f);
        mRocPlot.showCurrentGraph();
        final TextPoint2D data = mDataBundle.getMaxRangedPoint();
        if (data != null) {
          mStatusBar.setString(RocPlot.getMetricString(data.getY(), data.getX(), mDataBundle.getTotalVariants()) + " Threshold=" + data.getText());
        }
      }
    });

    final GridBagConstraints rangeConstraints = new GridBagConstraints();
    rangeConstraints.gridx = 0; rangeConstraints.gridy = 1;
    rangeConstraints.gridwidth = 2; rangeConstraints.fill = GridBagConstraints.HORIZONTAL;
    add(rangeSlider, rangeConstraints);
    setMaximumSize(new Dimension(getMaximumSize().width, getPreferredSize().height));
  }

  private static class CheckBoxListener implements ItemListener {
    private final RocLinePanel mPanel;

    public CheckBoxListener(RocLinePanel panel) {
      mPanel = panel;
    }
    @Override
    public void itemStateChanged(ItemEvent e) {
      mPanel.setSelected(mPanel.mCheckBox.isSelected());
    }
  }

  private String getTooltip(String path, String name) {
    return name + "=" + path;
  }

  public void setSelected(boolean selected) {
    mDataBundle.show(selected);
    mCheckBox.setSelected(selected);
    mRocPlot.showCurrentGraph();
  }

  /*public boolean isSelected() {
    return mData.get(mPath).show();
    //return mCheckBox.isSelected();
  }*/

  public String getText() {
    return mPath;
  }

  public void addActionListener(ActionListener listener) {
    mUpButton.addActionListener(listener);
    mDownButton.addActionListener(listener);
  }
}
