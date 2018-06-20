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
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.TreeSet;

import javax.imageio.ImageIO;

import com.rtg.util.HtmlReportHelper;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Generates report for Species module.
 *
 * Assumptions:
 * All reportable values are contiguous from column 0, and parsable as longs.
 * Last column is the name.
 *
 */
public final class SpeciesReport implements Report {

  static {
    System.setProperty("java.awt.headless", "true");
  }

  /** file name of the species output */
  public static final String SPECIES_TSV = "species.tsv";

  private static class SpeciesLine {
    private final float[] mValues;
    private boolean mHasValues = false;

    SpeciesLine(String[] parts, int length) {
      mValues = new float[length];
      for (int i = 0; i < mValues.length; ++i) {
//        System.err.println(parts.length + " : " + mValues.length + " : " + i + " : " + parts[i]);
        try {
          mValues[i] = Float.parseFloat(parts[i]);
        } catch (NumberFormatException e) {
          mValues[i] = 0;
        }
        if (mValues[i] >= 0.0005f) {
          mHasValues = true;
        }
      }
    }

    float getValue(int index) {
      return mValues[index];
    }

    boolean hasValues() {
      return mHasValues;
    }
  }

  private static class SpeciesData {
    final HashMap<String, SpeciesLine> mData = new HashMap<>();

    SpeciesData() { }

    private void put(String[] parts, boolean renamed, int colsToTrim) {
      mData.put(parts[parts.length - 1], new SpeciesLine(parts, parts.length - (renamed ? colsToTrim + 1 : colsToTrim)));
    }

    private float getValue(String species, int index) {
      final SpeciesLine sl = mData.get(species);
      return sl == null ? 0.0f : sl.getValue(index);
    }

    private boolean hasValues(String species) {
      final SpeciesLine sl = mData.get(species);
      return sl != null && sl.hasValues();
    }
  }

  @Override
   public void generateReport(ReportType type, File outputDir, File... inputDirs) throws IOException {
    if (type != ReportType.HTML) {
      throw new UnsupportedOperationException("Only " + ReportType.HTML + " reports implemented so far");
    }
    if (inputDirs == null || outputDir == null) {
      throw new NullPointerException();
    }

    final HtmlReportHelper hrh = new HtmlReportHelper(outputDir, "index");
    hrh.copyResources(ReportUtils.resourceArray("rtg.css", "table.js",
                      "table.css", "07_ascending.gif", "07_descending.gif")); //copy resources up here to create the report files sub dir as well

    final StringBuilder sampleNameMap = new StringBuilder();
    sampleNameMap.append("<table class=\"example\"><thead><tr><th>Sample</th><th>Location</th></tr></thead><tbody>");

    final int numSamples = inputDirs.length;
    final ArrayList<SpeciesData> data = new ArrayList<>();
    final ArrayList<String> columnNames = new ArrayList<>();
    final TreeSet<String> speciesNames = new TreeSet<>();
    final String[] sampleNames = new String[numSamples];

    Integer colsToTrim = 1;
    for (int sampleNum = 0; sampleNum < numSamples; ++sampleNum) {
      final File inputDir = inputDirs[sampleNum];
      sampleNames[sampleNum] = inputDirs[sampleNum].getName();
      sampleNameMap.append("<tr><td>").append(sampleNames[sampleNum]).append("</td><td>").append(inputDirs[sampleNum].getPath()).append("</td></tr>\n");

      // read data
      String line;
      boolean rename = true;
      File tsvFile = new File(inputDir, "species_rename.tsv");
      if (!tsvFile.exists()) {
        rename = false;
        tsvFile = new File(inputDir, SPECIES_TSV);
      }
      if (!tsvFile.exists()) {
        Diagnostic.warning("Could not find a species.tsv file in " + inputDir);
      } else {
        final SpeciesData sd = new SpeciesData();


        String[] lastCommentLine = null;
        boolean isFirstDataLine = true;
        try (BufferedReader br = new BufferedReader(new FileReader(tsvFile))) {
          while ((line = br.readLine()) != null) {
            final String[] parts = StringUtils.split(line, '\t');
            if (line.startsWith("#")) {
              if (line.startsWith("#Version")) {
                if (line.contains("species v2.0")) {
                  colsToTrim = 4;
                }
              } else {
                lastCommentLine = parts;
              }
            } else {
              if (isFirstDataLine) {  //get the column names from the last comment line (which is supposed to be the header)
                if (lastCommentLine == null) {
                  throw new IllegalArgumentException("Could not find header line in " + tsvFile);
                }
                if (sampleNum == 0) { //only grab the header from the first sample...
                  for (int headerCol = 0; headerCol < lastCommentLine.length; ++headerCol) {
                    columnNames.add(headerCol == 0 ? lastCommentLine[headerCol].substring(1) : lastCommentLine[headerCol]);
                  }
                  lastCommentLine = null;
                }
                isFirstDataLine = false;
              }
              if (parts.length != columnNames.size()) { //all tsv files must have same list of columns, apparently
                throw new IllegalArgumentException("Different number of columns in " + tsvFile + ": " + line);
              }
              final String name = parts[parts.length - 1];
              sd.put(parts, rename, colsToTrim);
              speciesNames.add(name);
            }
          }
        }
        data.add(sd);
      }

      if (columnNames.size() < 5) {
        throw new IOException("File does not appear to be a valid species.tsv file: " + tsvFile);
      }
    }
    sampleNameMap.append("</tbody></table>\n<br />\n");

    final TreeSet<String> populatedSpecies = new TreeSet<>();
    final float[] maxes = new float[columnNames.size() - colsToTrim];

    for (String speciesName : speciesNames) {
      boolean hasValues = false;
      for (int s = 0; s < numSamples; ++s) {
        final SpeciesData sd = data.get(s);
        if (sd != null && sd.hasValues(speciesName)) {
          hasValues = true;
        }
      }
      if (hasValues) {
        populatedSpecies.add(speciesName);

        for (int c = 0; c < columnNames.size() - colsToTrim; ++c) {
          for (int s = 0; s < numSamples; ++s) {
            final float value = data.get(s).getValue(speciesName, c);
            if (value > maxes[c]) {
              maxes[c] = value;
            }
          }
        }
      }
    }

    final String table = writeTable(data, columnNames, populatedSpecies, numSamples, colsToTrim, sampleNames, maxes, hrh);

    // write html
    final HashMap<String, String> replacements = new HashMap<>();
    replacements.put("__TITLE__", (numSamples > 1 ? "Multisample " : "") + "Species Report");
    replacements.put("__BODY_TEXT__", sampleNameMap.toString() + table);
    replacements.put("__RESOURCE_DIR__", hrh.getResourcesDirName() + "/");
    ReportUtils.writeHtml(ReportUtils.TEMPLATE_DIR + "/default.html", hrh.getReportFile(), replacements);
  }

  private String writeTable(List<SpeciesData> data, List<String> columnNames, TreeSet<String> populatedSpecies, int numSamples, int colsToTrim, String[] sampleNames, float[] maxes, HtmlReportHelper hrh) throws IOException {
    final StringBuilder table = new StringBuilder();
    table.append("<table class=\"example sort07 table-autosort:0\" id=\"page\" border=\"0\"><thead><tr>");
    for (int c = 0; c < columnNames.size() - colsToTrim; ++c) {
      final String colName = columnNames.get(c);
      table.append("<th class=\"table-sortable:numeric\" colspan=\"2\">").append(colName).append("</th>");
    }
    table.append("<th class=\"table-sortable:text\">Name</th>");
    table.append("</tr></thead>\n<tbody>");
    for (String speciesName : populatedSpecies) {
      table.append("<tr>");
      for (int c = 0; c < columnNames.size() - colsToTrim; ++c) {

        final float[] values = new float[numSamples];
        float max = -1.0f;
        final StringBuilder title = new StringBuilder();
        for (int s = 0; s < numSamples; ++s) {
          final SpeciesData sd = data.get(s);
          final float value = sd == null ? 0.0f : sd.getValue(speciesName, c);
          values[s] = value;
          if (value > max) {
            max = value;
          }
          if (title.length() != 0) {
            title.append("\n");
          }
          title.append(String.format(Locale.ROOT, "%1.6f", value)).append(" : ").append(sampleNames[s]);
        }
        table.append("<td>").append(String.format(Locale.ROOT, "%1.4f", max)).append("</td>");

        final String imageFile = createSparkBarImage(hrh.getResourcesDir(), maxes[c], values);
        table.append("<td><img src=\"").append(hrh.getResourcesDirName()).append('/').append(imageFile).append("\" title=\"").append(title.toString()).append("\"></td>");
      }
      table.append("<td>").append(speciesName).append("</td>");
      table.append("</tr>\n");
    }
    table.append("</tbody></table>");
    return table.toString();
  }

  private String floatToString(int[] values) {
    final StringBuilder sb = new StringBuilder();
    for (float v : values) {
      sb.append(v).append(" ");
    }
    return sb.toString();
  }

  private int mSparkBarCount = 0;
  private final HashMap<String, String> mSparkBarFiles = new HashMap<>();
  private static final int MAX_BAR_WIDTH = 64;
  private static final int GAP_WIDTH = 1;
  private static final int IMAGE_HEIGHT = 30;

  private String createSparkBarImage(File outputDir, float max, float[] values) throws IOException {
    float barMax = max;
    if (max <= 0.0f) {
      for (float f : values) {
        if (f > barMax) {
          barMax = f;
        }
      }
    }

    int barWidth = MAX_BAR_WIDTH / values.length - GAP_WIDTH;
    if (barWidth < 1) {
      barWidth = 1;
    }
    final int width = (barWidth + GAP_WIDTH) * values.length;
    final int height = IMAGE_HEIGHT;


    final int[] barValues = new int[values.length];
    for (int i = 0; i < values.length; ++i) {
      int val = (int) (height * values[i] / barMax);
      if (val > height) {
        val = height;
      }
      if (val < 1) {
        val = 1;
      }
      barValues[i] = val;
    }

    final String vstr = floatToString(barValues);
    if (mSparkBarFiles.containsKey(vstr)) {
      return mSparkBarFiles.get(vstr);
    }

    final BufferedImage bi = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
    final Graphics2D g = bi.createGraphics();
    g.setColor(Color.WHITE);
    g.fillRect(0, 0, width, height);

    g.setColor(Color.BLACK);
    g.fillRect(0, height - 1, width, 1);

    for (int i = 0; i < barValues.length; ++i) {
      final int h = barValues[i];
      g.setColor(new Color(255 - h * 255 / height, h * 255 / height, 0));
      g.fillRect((barWidth + GAP_WIDTH) * i, height - h, barWidth, h);
    }

    final String fileName = "sparkbar" + mSparkBarCount++ + ".png";
    final File output = new File(outputDir, fileName);
    try (FileOutputStream fos = new FileOutputStream(output)) {
      ImageIO.write(bi, "png", fos);
    }
    mSparkBarFiles.put(vstr, fileName);
    return fileName;
  }
}
