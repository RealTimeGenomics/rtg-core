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

package com.rtg.report;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Map.Entry;
import java.util.TreeMap;

import com.reeltwo.plot.Axis;
import com.reeltwo.plot.Box2D;
import com.reeltwo.plot.BoxPlot2D;
import com.reeltwo.plot.Datum2D;
import com.reeltwo.plot.Edge;
import com.reeltwo.plot.FillablePlot2D.FillStyle;
import com.reeltwo.plot.Graph2D;
import com.reeltwo.plot.KeyPosition;
import com.reeltwo.plot.LabelFormatter;
import com.reeltwo.plot.Plot2D;
import com.reeltwo.plot.Point2D;
import com.reeltwo.plot.PointPlot2D;
import com.reeltwo.plot.TextPlot2D;
import com.reeltwo.plot.TextPoint2D;
import com.reeltwo.plot.renderer.GraphicsRenderer;
import com.reeltwo.plot.ui.ImageWriter;
import com.rtg.launcher.CommonFlags;
import com.rtg.ngs.MapReportData;
import com.rtg.ngs.MapReportData.DistributionType;
import com.rtg.sam.RecordIterator;
import com.rtg.sam.SamFilterParams;
import com.rtg.sam.SamReadingContext;
import com.rtg.sam.SamRecordPopulator;
import com.rtg.sam.SamUtils;
import com.rtg.sam.ThreadedMultifileIterator;
import com.rtg.util.HistogramWithNegatives;
import com.rtg.util.HtmlReportHelper;
import com.rtg.util.License;
import com.rtg.util.MathUtils;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;

import htsjdk.samtools.SAMRecord;

/**
 */
public class MapReport extends MapSummaryReport {

  private static final int LINE_GRAPH_THRESHOLD = 10;

  private static class SamBamFilter implements FilenameFilter {
    @Override
    public boolean accept(File dir, String name) {
      final String caseInsensitiveName = name.toLowerCase(Locale.getDefault());
      return caseInsensitiveName.endsWith(".sam.gz") || caseInsensitiveName.endsWith(".sam") || caseInsensitiveName.endsWith(".bam");
    }
  }

  private final SamFilterParams mFilter;

  /**
   * Constructor for map report from existing data
   */
  public MapReport() {
    super();
    mFilter = null;
  }

  /**
   * Constructor for map report, allowing for report data regeneration
   * @param params the SAM filtering parameters (for other tools input SAM)
   */
  public MapReport(SamFilterParams params) {
    super();
    mFilter = params;
  }

  private void ensureReportData(File[] inputDirs) throws IOException {
    for (final File inputDir : inputDirs) {
      final File reportFile = new File(inputDir, MapReportData.MAP_REPORT_FILE_NAME);
      if (!(reportFile.exists() && checkVersion(reportFile))) {
        Diagnostic.userLog("Generating report data file " + reportFile + " from SAM files.");
        final MapReportData reportData = new MapReportData();
        final List<File> samFiles = Arrays.asList(FileUtils.listFiles(inputDir, new SamBamFilter()));
        if (samFiles.size() > 0) {
          if (mFilter == null) {
            throw new NoTalkbackSlimException("Report data needs to be (re)generated, but this has not been enabled");
          } else {
            final SamReadingContext context = new SamReadingContext(samFiles, 1, mFilter, SamUtils.getUberHeader(null, samFiles), null); // No need for a reference here unless we allow map to directly output CRAM
            try (final RecordIterator<SAMRecord> it = new ThreadedMultifileIterator<>(context, new SingletonPopulatorFactory<>(new SamRecordPopulator()))) {
              while (it.hasNext()) {
                final SAMRecord next = it.next();
                reportData.processRead(next);
              }
              if (it.getFilteredRecordsCount() + it.getInvalidRecordsCount() > 0) {
                Diagnostic.warning(it.getFilteredRecordsCount() + it.getInvalidRecordsCount() + " SAM records ignored due to filters");
              }
            }
          }
        }
        try (PrintStream ps = new PrintStream(reportFile)) {
          reportData.write(ps);
        }
      }
    }
  }

  @Override
  void reportContent(HtmlReportHelper helper, StringBuilder body, File[] inputDirs) throws IOException {

    ensureReportData(inputDirs);

    final MapReportData.Merger merger = new MapReportData.Merger();
    for (final File inputDir : inputDirs) {
      final File reportFile = new File(inputDir, MapReportData.MAP_REPORT_FILE_NAME);
      if (reportFile.exists() && checkVersion(reportFile)) {
        try (InputStream is = new FileInputStream(reportFile)) {
          merger.createMapReportData(is);
        }
      } else {
        throw new NoTalkbackSlimException("No compatible report data file: " + reportFile);
      }
    }
    final MapReportData reporter = merger.blendReportData();
    setCommandLines(reporter.getCommandLines());
    reporter.write(new File(helper.getBaseDir(), MapReportData.MAP_REPORT_FILE_NAME));

    super.reportContent(helper, body, inputDirs);

    if (inputDirs.length == 1) {
      body.append(reportSummary(helper, new File(inputDirs[0], CommonFlags.SUMMARY_FILE), 0));
    }
    if (inputDirs.length > 0) {
      if (reporter.isPairedEnd()) {
        if (License.isDeveloper()) {
          body.append(reportDistribution(helper, reporter, DistributionType.MAPC.longName(), DistributionType.MAPC));
        }
        for (DistributionType dt : MapReportData.getPairedBaseTypes()) {
            body.append(reportDistribution(helper, reporter, dt.longName() + " (Mated Only)", dt.getType(false, true), dt.getType(true, true)));
            body.append(reportDistribution(helper, reporter, dt.longName() + " (Unmated Only)", dt.getType(false, false), dt.getType(true, false)));
        }
        final DistributionType[] others = {DistributionType.MORI, DistributionType.FLEN, DistributionType.RLENU, DistributionType.MAPQ};
        for (DistributionType dt : others) {
            body.append(reportDistribution(helper, reporter, dt.longName(), dt));
        }
      } else {
        for (DistributionType dt : DistributionType.values()) {
            body.append(reportDistribution(helper, reporter, dt.longName(), dt));
        }
      }
      if (inputDirs.length > 1) {
        for (int i = 0; i < inputDirs.length; ++i) {
          body.append(reportSummary(helper, new File(inputDirs[i], CommonFlags.SUMMARY_FILE), i));
        }
      }
    }
  }

  /**
   * Check the version number of the report file is the same as the current MapReport version.
   * @param reportFile the report file to check the version of.
   * @return true if the version matches, false otherwise.
   * @throws IOException
   */
  private boolean checkVersion(File reportFile) throws IOException {
    try (InputStream is = new FileInputStream(reportFile)) {
      try (BufferedReader reader = new BufferedReader(new InputStreamReader(is))) {
        String line;
        while ((line = reader.readLine()) != null) {
          if (line.startsWith("#Version") && line.contains(MapReportData.VERSION)) {
            return true;
          }
        }
      }
    }
    return false;
  }

  private HistogramWithNegatives[] getHistogramsForTypes(MapReportData reporter, DistributionType... dt) {
    final HistogramWithNegatives[] hists = new HistogramWithNegatives[dt.length];
    for (int i = 0; i < dt.length; ++i) {
      hists[i] = reporter.getHistogram(dt[i]);
    }
    return hists;
  }

  private String reportSummary(HtmlReportHelper helper, File summaryFile, int fileIndex) throws IOException {
    if (!summaryFile.exists()) {
      return "";
    }
    final boolean isPaired = FileUtils.fileToString(summaryFile).contains("mated");
    final HistogramWithNegatives[] grams = new HistogramWithNegatives[isPaired ? 2 : 1];
    for (int i = 0; i < grams.length; ++i) {
      grams[i] = new HistogramWithNegatives();
    }
    final StringBuilder summaryText = new StringBuilder();
    //summaryText.append("<html><body><pre>\n");
    final ArrayList<String> titles = new ArrayList<>();
    try (BufferedReader reader = new BufferedReader(new FileReader(summaryFile))) {
      String line;
      final int numParts = isPaired ? 5 : 3;
      while ((line = reader.readLine()) != null) {
        final String trimmedLine = line.trim();
        if (trimmedLine.length() > 0) {
          summaryText.append(line).append("\n");
          final String[] parts = trimmedLine.split("\\s+", numParts);
          if (parts.length == numParts) {
            final String title = parts[numParts - 1];
            if (title.startsWith("ma") || title.startsWith("un")) {
              for (int i = 0; i < grams.length; ++i) {
                grams[i].increment(titles.size(), Long.parseLong(parts[i]));
              }
              final String[] tParts = title.split("\\s+");
              if (tParts[tParts.length - 1].endsWith(")")) {
                titles.add("  " + tParts[0] + " " + tParts[tParts.length - 3] + " " + tParts[tParts.length - 2] + " " + tParts[tParts.length - 1]);
              } else {
                titles.add("  " + tParts[0] + " " + tParts[tParts.length - 2] + " " + tParts[tParts.length - 1]);
              }
            }
          }
        }
      }
    }
    //summaryText.append("</pre></body></html>\n");

    final DataTable dTable = new DataTable(titles.toArray(new String[titles.size()]), grams);
    final String sectionTitle = "Summary";
    final File png = new File(helper.getResourcesDir(), sectionTitle + "_" + fileIndex +  ".png");
    dTable.toGraphImage(png, PlotType.BOX);
    //final File sOutFile = FileUtils.stringToFile(summaryText.toString(), new File(outputDir, sectionTitle + "_" + fileIndex + ".html"));
    //return getHtmlChunk(sectionTitle + " for " + summaryFile.getParent(), png, sOutFile);
    return getHtmlChunk(helper, sectionTitle + " for " + summaryFile.getParent(), png, summaryText.toString());
  }

  private String reportDistribution(HtmlReportHelper helper, MapReportData reporter, String sectionTitle, DistributionType... dt) throws IOException {
    if (!dt[0].showData() && !dt[0].showImage()) { //no report to be shown for this type at this time
      return "";
    }
    final HistogramWithNegatives[] grams = getHistogramsForTypes(reporter, dt);
    boolean empty = true;
    for (HistogramWithNegatives gram : grams) {
      if (gram.getLength() > 0) {
        empty = false;
        break;
      }
    }
    if (empty) {
      return "";
    }
    final DataTable dTable = new DataTable(grams);
    final PlotType pType;
    if (dt[0].equals(DistributionType.MORI)) {
      pType = PlotType.MORI;
    } else if (dt[0].ordinal() >= DistributionType.ORI.ordinal() && dt[0].ordinal() <= DistributionType.ORI2M.ordinal()) {
      pType = PlotType.ORI;
    } else if (dt[0].equals(DistributionType.MAPC)) {
      pType = PlotType.MAPC;
    } else if (dt[0].equals(DistributionType.MAPQ)) {
      pType = PlotType.BOX;
    } else if (dt[0].equals(DistributionType.FLEN) || dTable.xWidth() > LINE_GRAPH_THRESHOLD) {
      pType = PlotType.LINE;
    } else {
      pType = PlotType.BOX;
    }
    File png = null;
    String table = null;
    if (dt[0].showImage()) {
      png = new File(helper.getResourcesDir(), sectionTitle.replace(' ', '-') + ".png");
      dTable.toGraphImage(png, pType);
    }
    if (dt[0].showData()) {
      table = dTable.toHtmlTable(pType);
    }
    return getHtmlChunkF(helper, sectionTitle, png, table);
  }

  private String getHtmlChunkF(HtmlReportHelper helper, final String sectionTitle, final File image, final String text) {
    final StringBuilder html = new StringBuilder();
    html.append("<h3 style=\"clear: both\">").append(sectionTitle).append("</h3>").append(StringUtils.LS);
    html.append("<div style=\"max-width: 954px; overflow: auto\">").append(StringUtils.LS);
    if (image != null) {
      html.append("<div style=\"float: left; padding-right: 5px\"><img src=\"").append(helper.getResourcesDirName()).append("/").append(image.getName()).append("\"></div>").append(StringUtils.LS);
    }
    if (text != null) {
      html.append(text).append(StringUtils.LS);
    }
    html.append("</div>").append(StringUtils.LS);
    return html.toString();
  }

  private String getHtmlChunk(HtmlReportHelper helper, final String sectionTitle, final File image, final String text) {
    return "<h3>" + sectionTitle + "</h3>" + StringUtils.LS + "<table><tr>" + StringUtils.LS + "<td valign=\"top\"><img src=\"" + helper.getResourcesDirName() + "/" + image.getName() + "\"></td>" + StringUtils.LS + "<td valign=\"top\"><pre>" + helper.getResourcesDirName() + "/" + text + "</pre></td>" + StringUtils.LS + "</tr></table>" + StringUtils.LS;
  }

  /**
   * Types of plot available
   */
  protected enum PlotType {
    BOX,
    LINE,
    ORI,
    MORI,
    MAPC
  }

  private static class StringLabelFormatter implements LabelFormatter {
    private final String[] mLabels;

    StringLabelFormatter(String[] labels) {
      mLabels = labels;
    }

    @Override
    public String format(float f) {
      if (Float.compare((float) MathUtils.round(f), f) != 0) {
        return "";
      }
      return mLabels[(int) f % mLabels.length];
    }
  }

  private static class DataTable {

    private static class Values {
      private final long[] mValues;
      Values(int size) {
        mValues = new long[size];
      }
      void setValue(int index, long value) {
        mValues[index] = value;
      }
      long getValue(int index) {
        return mValues[index];
      }
    }

    final TreeMap<Integer, Values> mMap = new TreeMap<>();
    final int mMin;
    final int mMax;
    final int mSize;
    static final float WIDTH = 0.8f;
    static final float H_WIDTH = WIDTH / 2.0f;
    final float mBarSize;
    static final int[] COLORS = {2, 1};
    static final String[] KEY_TITLES = {"left", "right"};
    static final String[] ORI_TITLES = {"forward", "reverse"};
    static final String[] MORI_TITLES = {"FF", "FR", "RF", "RR"};
    static final String[] MAPPING_TITLES = {"Mapped", "Mated", "Unmated", "Unmapped"};
    final String[] mTitles;

    DataTable(String[] titles, HistogramWithNegatives... grams) {
      assert grams.length > 0 && grams.length <= 2;
      int min = Integer.MAX_VALUE;
      int max = 0;
      for (int j = 0; j < grams.length; ++j) {
        final HistogramWithNegatives gram = grams[j];
        for (int i = gram.min(); i < gram.max(); ++i) {
          final long val = gram.getValue(i);
          if (val != 0) {
            Values vals = mMap.get(i);
            if (vals == null) {
              vals = new Values(grams.length);
              mMap.put(i, vals);
            }
            vals.setValue(j, val);
            max = Math.max(max, i);
            min = Math.min(min, i);
          }
        }
      }
      mMin = min;
      mMax = max;
      mSize = grams.length;
      mBarSize = WIDTH / mSize;
      mTitles = titles;
    }

    DataTable(HistogramWithNegatives... grams) {
      this(null, grams);
    }

    int xWidth() {
      return mMax - mMin;
    }

    LabelFormatter getLabelFormatter(PlotType type) {
      switch (type) {
        case MORI:
          return new StringLabelFormatter(MORI_TITLES);
        case ORI:
          return new StringLabelFormatter(ORI_TITLES);
        case MAPC:
          return new StringLabelFormatter(MAPPING_TITLES);
        default:
          return null;
      }
    }

    String toHtmlTable(PlotType type) {
      final LabelFormatter formatter = getLabelFormatter(type);
      final StringBuilder sb = new StringBuilder();
      sb.append("<div style=\"overflow:auto; max-height:240px;max-width:310px;border: 1px solid black\"><table>").append(StringUtils.LS);
      for (Entry<Integer, Values> e : mMap.entrySet()) {
        sb.append("<tr><td align=\"right\">");
        if (formatter != null) {
          sb.append(formatter.format(e.getKey()));
        } else {
          sb.append(e.getKey());
        }
        sb.append("</td>");
        for (int i = 0; i < mSize; ++i) {
          sb.append("<td align=\"right\">");
          sb.append(e.getValue().getValue(i));
          sb.append("</td>");
        }
        sb.append("</tr>").append(StringUtils.LS);
      }
      sb.append("</table></div>").append(StringUtils.LS);
      return sb.toString();
    }

    Plot2D getPlot(PlotType type) {
      switch (type) {
        case LINE:
          final PointPlot2D pPlot = new PointPlot2D();
          pPlot.setLines(true);
          pPlot.setPoints(false);
          pPlot.setLineWidth(2);
          return pPlot;
        default:
          final BoxPlot2D bPlot = new BoxPlot2D();
          bPlot.setFill(FillStyle.COLOR);
          return bPlot;
      }
    }


    Datum2D getDatum(PlotType type, int plotIndex, int x, long y) {
      final Datum2D datum;
      switch (type) {
        case LINE:
          datum = new Point2D(x, y);
          break;
        default:
          datum = new Box2D(x - H_WIDTH + mBarSize * plotIndex, 0.0f, x - H_WIDTH + mBarSize * (plotIndex + 1), y);
          break;
      }
      return datum;
    }

    Datum2D[] getDatumArray(PlotType type, int size) {
      switch (type) {
        case LINE:
          return new Point2D[size];
        default:
          return new Box2D[size];
      }
    }

    void toGraphImage(File output, PlotType type) throws IOException {
      final Plot2D[] plot = new Plot2D[mSize];
      final Datum2D[][] data = new Datum2D[mSize][];
      for (int i = 0; i < mSize; ++i) {
        plot[i] = getPlot(type);
        plot[i].setColor(COLORS[i]);
        plot[i].setTitle(KEY_TITLES[i]);
        data[i] = getDatumArray(type, mMap.size());
      }
      int j = 0;
      for (Entry<Integer, Values> e : mMap.entrySet()) {
        for (int i = 0; i < mSize; ++i) {
          data[i][j] = getDatum(type, i, e.getKey(), e.getValue().getValue(i));
        }
        ++j;
      }

      final Graph2D graph = new Graph2D();
      switch (type) {
        case MORI:
        case MAPC:
          graph.setRange(Axis.X, -0.5f, 3.5f);
          break;
        case ORI:
          graph.setRange(Axis.X, -0.5f, 1.5f);
          break;
        default:
          graph.setRange(Axis.X, mMin - 0.5f, mMax + 0.5f);
          break;
      }
      if (xWidth() <= LINE_GRAPH_THRESHOLD) {
        graph.setTic(Graph2D.X, Graph2D.ONE, 1.0f);
      }
      final LabelFormatter lf = getLabelFormatter(type);
      if (lf != null) {
        graph.setTicLabelFormatter(Axis.X, Edge.MAIN, lf);
      }
      graph.setGrid(true);
      if (mSize > 1) {
        graph.setShowKey(true);
        graph.setKeyHorizontalPosition(KeyPosition.OUTSIDE);
        graph.setKeyVerticalPosition(KeyPosition.TOP);
      } else {
        graph.setShowKey(false);
      }
      for (int i = 0; i < mSize; ++i) {
        plot[i].setData(data[i]);
        graph.addPlot(plot[i]);
      }
      if (mTitles != null) {
        final TextPlot2D tPlot = new TextPlot2D();
        final TextPoint2D[] titles = new TextPoint2D[mTitles.length];
        for (int i = 0; i < mTitles.length; ++i) {
          titles[i] = new TextPoint2D(i, 0, mTitles[i]);
        }
        tPlot.setData(titles);
        tPlot.setVertical(true);
        graph.addPlot(tPlot);
        graph.setGrid(Axis.X, Edge.MAIN, false);
        graph.setTicLabelFormatter(Axis.X, Edge.MAIN, new StringLabelFormatter(new String[] {""}));
      }

      final ImageWriter iw = new ImageWriter(new GraphicsRenderer(new Color[] {Color.RED, new Color(0, 192, 0), new Color(0, 128, 255)}));
      iw.toPNG(output, graph, 640, 240, null);
    }
  }
}
