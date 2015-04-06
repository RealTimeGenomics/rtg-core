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
import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.File;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.CountDownLatch;

import javax.swing.AbstractAction;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JSplitPane;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import com.reeltwo.jumble.annotations.JumbleIgnore;
import com.reeltwo.plot.Graph2D;
import com.reeltwo.plot.KeyPosition;
import com.reeltwo.plot.Plot2D;
import com.reeltwo.plot.Point2D;
import com.reeltwo.plot.PointPlot2D;
import com.reeltwo.plot.renderer.Mapping;
import com.reeltwo.plot.ui.PlotPanel;
import com.reeltwo.plot.ui.ZoomPlotPanel;
import com.rtg.util.Resources;
import com.rtg.util.io.FileUtils;

/**
 * Starts a new Swing window for displaying {@code Graph2D}s in. The window has
 * zooming and picture in picture functionality enabled.
 *
 */
@JumbleIgnore
public final class RocPlot {


  /** Minimum allowed line width */
  public static final int LINE_WIDTH_MIN = 1;
  /** Maximum allowed line width */
  public static final int LINE_WIDTH_MAX = 10;
  
  private final ProgressBarDelegate mProgressBarDelegate;

  private final JPanel mMainPanel;
  /** panel showing plot */
  private final PlotPanel mPlotPanel;
  private final RocZoomPlotPanel mZoomPP;
  /** a progress bar */
  private final JProgressBar mProgressBar;
  /** pop up menu */
  private final JPopupMenu mPopup;

  private final JLabel mIconLabel;

  private final RocLinesPanel mRocLinesPanel;
  private final JSlider mLineWidthSlider;
  private final JCheckBox mScoreCB;
  private final JCheckBox mSelectAllCB;
  private final JButton mOpenButton;
  private final JTextField mTitleEntry;
  private JSplitPane mSplitPane;
  private final JLabel mStatusLabel;

  // Graph data and state
  final Map<String, DataBundle> mData = Collections.synchronizedMap(new HashMap<String, DataBundle>());
  boolean mShowScores = true;
  int mLineWidth = 2;

  private float mMaxXHi = -1.0f;
  private float mMaxYHi = -1.0f;

  private final JScrollPane mScrollPane;


  /** Creates a new swing plot. */
  RocPlot() {
    mMainPanel = new JPanel();
    mPlotPanel = new PlotPanel(true);
    mZoomPP = new RocZoomPlotPanel(mPlotPanel, mMainPanel);
    mZoomPP.setOriginIsMin(true);
    mProgressBar = new JProgressBar(-1, -1);
    mProgressBarDelegate = new ProgressBarDelegate(mProgressBar);
    mStatusLabel = new JLabel();
    mPopup = new JPopupMenu();
    mRocLinesPanel = new RocLinesPanel(this);
    mScrollPane = new JScrollPane(mRocLinesPanel, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
    mScrollPane.setWheelScrollingEnabled(true);
    //System.err.println("scroll: " + mScrollPane.getPreferredSize());
 //   mRocLinesPanel.setLayout(new GridLayout(0, 1));
    mLineWidthSlider = new JSlider(JSlider.HORIZONTAL, LINE_WIDTH_MIN, LINE_WIDTH_MAX, 1);
    mScoreCB = new JCheckBox("Show Scores");
    mScoreCB.setSelected(true);
    mSelectAllCB = new JCheckBox("Select / Deselect all");
    mTitleEntry = new JTextField("ROC");
    mOpenButton = new JButton("+");
    final ImageIcon icon = createImageIcon("com/rtg/graph/resources/realtimegenomics_logo.png", "RTG Logo");
    mIconLabel = new JLabel(icon);
    mIconLabel.setBackground(new Color(16, 159, 205));
    mIconLabel.setForeground(Color.WHITE);
    mIconLabel.setOpaque(true);
    mIconLabel.setFont(new Font("Arial", Font.BOLD, 24));
    mIconLabel.setHorizontalAlignment(JLabel.LEFT);
    mIconLabel.setIconTextGap(50);
    configureUI();
  }

  protected static ImageIcon createImageIcon(String path, String description) {
    final java.net.URL imgURL = Resources.getResource(path);
    if (imgURL != null) {
      return new ImageIcon(imgURL, description);
    } else {
      System.err.println("Couldn't find file: " + path);
      return null;
    }
  }

  /**
   * Layout and show the GUI.
   */
  private void configureUI() {
    mMainPanel.setLayout(new BorderLayout());
    final JPanel pane = new JPanel(new BorderLayout());
    pane.add(mPlotPanel, BorderLayout.CENTER);
    final Box rightPanel = new Box(BoxLayout.Y_AXIS);
    mSplitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, pane, rightPanel);
    mSplitPane.setContinuousLayout(true);
    mSplitPane.setOneTouchExpandable(true);
    mSplitPane.setResizeWeight(1);

    mMainPanel.add(mSplitPane, BorderLayout.CENTER);
    mMainPanel.add(mStatusLabel, BorderLayout.SOUTH);

    mPopup.setLightWeightPopupEnabled(false);
    mPopup.add(mZoomPP.getZoomOutAction());

    mPopup.addSeparator();
    mPopup.add(mPlotPanel.getPrintAction());
    mPopup.add(mPlotPanel.getSaveImageAction());
    mPopup.add(mPlotPanel.getSnapShotAction());

    // listener to show popup
    mPlotPanel.addMouseListener(new PopupListener());

    mPopup.addSeparator();

    //Set up the content pane.
    mPlotPanel.setBackground(Color.WHITE);
    mPlotPanel.setGraphBGColor(new Color(0.8f, 0.9f, 1.0f), Color.WHITE);
    mPlotPanel.setGraphShadowWidth(4);

    final JPanel controlPanel = new JPanel(new GridLayout(0, 1));
    controlPanel.add(new JLabel("Line Width", JLabel.CENTER));
    mLineWidthSlider.setSnapToTicks(true);
    mLineWidthSlider.setValue(mLineWidth);
    //mLineWidthSlider.setPreferredSize(new Dimension(300, mLineWidthSlider.getPreferredSize().height));
    mLineWidthSlider.addChangeListener(new ChangeListener() {
      @Override
      public void stateChanged(ChangeEvent e) {
        mLineWidth = mLineWidthSlider.getValue();
        showCurrentGraph();
      }
    });
    controlPanel.add(mLineWidthSlider);

    mScoreCB.addItemListener(new ItemListener() {
      @Override
      public void itemStateChanged(ItemEvent e) {
        mShowScores = mScoreCB.isSelected();
        showCurrentGraph();
      }
    });
    controlPanel.add(mScoreCB);

    final JPanel titlePanel = new JPanel(new GridLayout(0, 1));
    titlePanel.add(new JLabel("Title", JLabel.CENTER));
    mTitleEntry.addKeyListener(new KeyAdapter() {
      @Override
      public void keyPressed(KeyEvent e) {
        if (e.getKeyCode() == KeyEvent.VK_ENTER) {
          mIconLabel.setText(mTitleEntry.getText());
          showCurrentGraph();
        }
      }

    });
    titlePanel.add(mTitleEntry);
    controlPanel.add(titlePanel);

    rightPanel.add(controlPanel);

    final JPanel checkControlPanel = new JPanel(new GridBagLayout());
    final GridBagConstraints openButtonConstraints = new GridBagConstraints();
    openButtonConstraints.gridx = 0; openButtonConstraints.gridy = 0;
    openButtonConstraints.anchor = GridBagConstraints.LINE_START;
    checkControlPanel.add(mOpenButton, openButtonConstraints);
    final GridBagConstraints selectAllConstraints = new GridBagConstraints();
    selectAllConstraints.gridx = 0; selectAllConstraints.gridy = 1;
//    selectAllConstraints.gridx = 2;
    checkControlPanel.add(mSelectAllCB, selectAllConstraints);

    mSelectAllCB.addItemListener(new ItemListener() {
      @Override
      public void itemStateChanged(ItemEvent e) {
        for (final Component component : mRocLinesPanel.getComponents()) {
          final RocLinePanel cp = (RocLinePanel) component;
          cp.setSelected(mSelectAllCB.isSelected());
        }
      }
    });
    mSelectAllCB.setSelected(true);
    mOpenButton.addActionListener(new ActionListener() {
      private File mFileChooserParent = null;

      @Override
      public void actionPerformed(ActionEvent e) {
        final JFileChooser chooser = new JFileChooser();
        if (mFileChooserParent != null) {
          chooser.setCurrentDirectory(mFileChooserParent);
        }
        chooser.showOpenDialog(checkControlPanel);
        final File f = chooser.getSelectedFile();
        if (f != null) {
          try {
            mFileChooserParent = f.getParentFile();
            loadFile(f, f.getCanonicalFile().getParentFile().getName(), true);
            updateProgress();
            showCurrentGraph();
          } catch (final IOException e1) {
            JOptionPane.showInternalMessageDialog(checkControlPanel, "Could not open file " + f.getPath());
          }
        }
      }
    });
    final JPanel namePanel = new JPanel(new GridBagLayout());
    final GridBagConstraints checkConstraints = new GridBagConstraints();
    checkConstraints.gridx = 0; checkConstraints.gridy = 0;
    checkConstraints.anchor = GridBagConstraints.FIRST_LINE_START;
    namePanel.add(checkControlPanel, checkConstraints);
    final GridBagConstraints scrollConstraints = new GridBagConstraints();
    scrollConstraints.gridx = 0; scrollConstraints.gridy = 1;
    scrollConstraints.weightx = 2; scrollConstraints.weighty = 2;
    scrollConstraints.fill = GridBagConstraints.BOTH;
    mRocLinesPanel.setPreferredSize(new Dimension(mScrollPane.getViewport().getViewSize().width, mRocLinesPanel.getPreferredSize().height));
    namePanel.add(mScrollPane, scrollConstraints);
    rightPanel.add(namePanel);
//    System.err.println("Scroll: " + mScrollPane.getSize());

    pane.add(mProgressBar, BorderLayout.SOUTH);

    mIconLabel.setText(mTitleEntry.getText());
    pane.add(mIconLabel, BorderLayout.NORTH);
  }

  // Adds the notion of painting a current crosshair position
  @JumbleIgnore
  private class RocZoomPlotPanel extends ZoomPlotPanel {
    private final PlotPanel mPlotPanel;
    private Point mCrosshair; // In TP / FP coordinates.
    RocZoomPlotPanel(PlotPanel plotPanel, Container container) {
      super(plotPanel, container);
      mPlotPanel = plotPanel;
    }
    @Override
    public void paint(Graphics g) {
      super.paint(g);
      final Mapping[] mapping = mPlotPanel.getMapping();
      if (mapping != null && mapping.length > 1 && mCrosshair != null) {
        Point p = new Point((int) mapping[0].worldToScreen(mCrosshair.x), (int) mapping[1].worldToScreen(mCrosshair.y));
        p = SwingUtilities.convertPoint(mPlotPanel, p, this);
        g.setColor(Color.BLACK);
        final int size = 9;
        g.drawLine(p.x - size, p.y - size, p.x + size, p.y + size);
        g.drawLine(p.x - size, p.y + size, p.x + size, p.y - size);
      }
    }
    void setCrossHair(Point p) {
      mCrosshair = p;
    }
  }

  // Adds the notion of the baseline total number of variants, for calculating sensitivity
  @JumbleIgnore
  static class RocGraph2D extends Graph2D {
    private final int mMaxVariants;
    RocGraph2D(ArrayList<String> lineOrdering, int lineWidth, boolean showScores, Map<String, DataBundle> data, String title) {
      setKeyVerticalPosition(KeyPosition.BOTTOM);
      setKeyHorizontalPosition(KeyPosition.RIGHT);
      setGrid(true);
      setLabel(Graph2D.Y, "True Positives");
      setLabel(Graph2D.X, "False Positives");
      setTitle(title);

      int maxVariants = -1;
      for (int i = 0; i < lineOrdering.size(); i++) {
        final DataBundle db = data.get(lineOrdering.get(i));
        if (db.show()) {
          addPlot(db.getPlot(lineWidth, i));
          if (showScores) {
            addPlot(db.getScorePoints(lineWidth, i));
            addPlot(db.getScoreLabels());
          }
          if (db.getTotalVariants() > maxVariants) {
            maxVariants = db.getTotalVariants();
          }
        }
      }

      if (maxVariants > 0) {
        setRange(Graph2D.Y, 0, maxVariants);
        setTitle(title + " (baseline total = " + maxVariants + ")");

        setRange(Graph2D.Y, Graph2D.TWO, 0, 100);
        setShowTics(Graph2D.Y, Graph2D.TWO, true);
        setGrid(Graph2D.Y, Graph2D.TWO, false);
        setLabel(Graph2D.Y, Graph2D.TWO, "%");
        // dummy plot to show Y2 axis
        final PointPlot2D pp = new PointPlot2D(Graph2D.ONE, Graph2D.TWO);
        addPlot(pp);
      }
      mMaxVariants = maxVariants;
    }

    public int getMaxVariants() {
      return mMaxVariants;
    }
  }

  void showCurrentGraph() {
    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        final Graph2D graph = new RocGraph2D(RocPlot.this.mRocLinesPanel.plotOrder(), RocPlot.this.mLineWidth, RocPlot.this.mShowScores, RocPlot.this.mData, RocPlot.this.mTitleEntry.getText());
        maintainZoomMax(graph);
        graph.addPlot(invisibleGraph());
        mZoomPP.setGraph(graph, true);
      }
    });
  }

  private Plot2D invisibleGraph() {
    // Invisible graph to maintain graph size when no lines are shown
    final PointPlot2D plot = new PointPlot2D();
    plot.setData(Arrays.asList(new Point2D(0, 0), new Point2D(mMaxXHi, mMaxYHi)));
    plot.setLines(false);
    plot.setPoints(false);
    return plot;
  }

  private void maintainZoomMax(Graph2D graph) {
    mMaxXHi = Math.max(mMaxXHi, graph.getHi(Graph2D.X, Graph2D.ONE));
    mMaxYHi = Math.max(mMaxYHi, graph.getHi(Graph2D.Y, Graph2D.ONE));
  }

  void updateProgress() {
    mProgressBarDelegate.done();
  }

  public ProgressBarDelegate getProgressBarDelegate() {
    return mProgressBarDelegate;
  }


  /**
   * Set the title of the plot
   * @param title plot title
   */
  public void setTitle(final String title) {
    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        mTitleEntry.setText(title);
      }
    });
  }

  /**
   * Set whether to show scores on the plot lines
   * @param flag show scores
   */
  public void showScores(boolean flag) {
    mShowScores = flag;
    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        mScoreCB.setSelected(mShowScores);
      }
    });
  }

  /**
   * Set whether to show the open file button
   * @param flag show open file button
   */
  public void showOpenButton(boolean flag) {
    mOpenButton.setVisible(flag);
  }

  /**
   * Set the line width slider to the given value
   * @param width line width
   */
  public void setLineWidth(int width) {
    mLineWidth = width < LINE_WIDTH_MIN ? LINE_WIDTH_MIN : width > LINE_WIDTH_MAX ? LINE_WIDTH_MAX : width;
    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        mLineWidthSlider.setValue(mLineWidth);
      }
    });
  }

  /**
   * Sets the split pane divider location
   * @param loc proportional location
   */
  public void setSplitPaneDividerLocation(double loc) {
    mSplitPane.setDividerLocation(loc);
  }

  /**
   * Returns the main application panel
   * @return main panel
   */
  public JPanel getMainPanel() {
    return mMainPanel;
  }

  /**
   * Returns the zooming part of the plot panel
   * @return zoom plot panel
   */
  public ZoomPlotPanel getZoomPlotPanel() {
    return mZoomPP;
  }

  /**
   * Set a status message
   * @param message test to display
   */
  public void setStatus(String message) {
    mStatusLabel.setText(message);
  }

  private void loadData(ArrayList<File> files, ArrayList<String> names, boolean showProgress) throws IOException {
    for (int i = 0; i < files.size(); i++) {
      final File f = files.get(i);
      final String name = names.get(i);
      loadFile(f, name, showProgress);
    }
    if (showProgress) {
      updateProgress();
    }
  }

  private void loadFile(final File f, final String name, boolean showProgress) throws IOException {
    final DataBundle data = ParseRocFile.loadStream(mProgressBarDelegate, FileUtils.createInputStream(f, false), name, showProgress);
    addLine(f.getAbsolutePath(), data);
  }

  private void addLine(String path, DataBundle dataBundle) {
    mData.put(path, dataBundle);
    mRocLinesPanel.addLine(new RocLinePanel(this, path, dataBundle.getTitle(), dataBundle, mProgressBar));
    showCurrentGraph();
  }

  /**
   * A class required to listen for right-clicks
   */
  private class PopupListener extends MouseAdapter {
    @Override
    public void mouseClicked(MouseEvent e) {
      final Point p = e.getPoint();
      final Mapping[] mapping = mPlotPanel.getMapping();
      final int maxVariants = ((RocGraph2D) mPlotPanel.getGraph()).getMaxVariants();
      if (mapping != null && mapping.length > 1 && p != null) {
        final float x = mapping[0].screenToWorld((float) p.getX());
        final float y = mapping[1].screenToWorld((float) p.getY());
        if (x >= 0 && y >= 0 && (x + y > 0)) {
          String message = String.format("TP=%.0f FP=%.0f Precision=%.2f%%", y, x, y / (x + y) * 100);
          if (maxVariants > 0) {
            message += String.format(" Sensitivity=%.2f%%", y / maxVariants * 100);
          }
          mZoomPP.setCrossHair(new Point((int) x, (int) y));
          mProgressBar.setString(message);
        } else {
          mZoomPP.setCrossHair(null);
          mProgressBar.setString("");
        }
      }
    }

    @Override
    public void mousePressed(MouseEvent e) {
      maybeShowPopup(e);
    }

    @Override
    public void mouseReleased(MouseEvent e) {
      maybeShowPopup(e);
    }

    private void maybeShowPopup(MouseEvent e) {
      if (e.isPopupTrigger()) {
        mPopup.show(e.getComponent(), e.getX(), e.getY());
      }
    }
  }


  void rocStandalone(ArrayList<File> fileList, ArrayList<String> nameList, String title, boolean scores, final boolean hideSidePanel, int lineWidth) throws InterruptedException, InvocationTargetException, IOException {
    final RocPlot rp = new RocPlot();
    rp.setLineWidth(lineWidth);
    final JFrame frame = new JFrame("ROC");
    frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
    frame.setLayout(new BorderLayout());
    frame.add(rp.mMainPanel, BorderLayout.CENTER);
    frame.setGlassPane(rp.mZoomPP);
    frame.getGlassPane().setVisible(true);
    final CountDownLatch lock = new CountDownLatch(1);
    rp.mPopup.add(new AbstractAction("Exit", null) {
      @Override
      public void actionPerformed(ActionEvent e) {
        frame.setVisible(false);
        frame.dispose();
        lock.countDown();
      }
    });
    rp.showScores(scores);
    if (title != null) {
      rp.setTitle(title);
    }
    SwingUtilities.invokeAndWait(new Runnable() {
      @Override
      public void run() {
        frame.pack();
        frame.setSize(1024, 768);
        frame.setLocation(50, 50);
        frame.setVisible(true);
        rp.showCurrentGraph();
        if (hideSidePanel) {
          rp.setSplitPaneDividerLocation(1.0);
        }
      }
    });
    rp.loadData(fileList, nameList, true);
    SwingUtilities.invokeAndWait(new Runnable() {
      @Override
      public void run() {
        rp.mZoomPP.setGraph(new RocGraph2D(rp.mRocLinesPanel.plotOrder(), rp.mLineWidth, rp.mShowScores, rp.mData, rp.mTitleEntry.getText()), false);
      }
    });
    lock.await();
  }
}
