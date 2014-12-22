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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Locale;

import com.rtg.similarity.SimilaritySvd;
import com.rtg.util.HtmlReportHelper;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 */
public class SimilarityReport implements Report {

  private static class Heats {
    private final String mName;
    private final int[] mValues;

    Heats(String name, int count) {
      if (name == null) {
        throw new NullPointerException("Name cannot be null");
      }
      mName = name;
      mValues = new int[count];
    }

    String getName() {
      return mName;
    }

    void setValue(int index, int value) {
      mValues[index] = value;
    }

    int getValue(int index) {
      return mValues[index];
    }
  }


  @Override
  public void generateReport(ReportType type, File[] inputDirs, File outputDir) throws IOException {
    if (type != ReportType.HTML) {
      throw new UnsupportedOperationException("Only " + ReportType.HTML + " reports implemented so far");
    }
    if (inputDirs == null) {
      throw new NullPointerException();
    }
    if (inputDirs.length == 0) {
      throw new IllegalArgumentException("No input directories given.");
    }
    if (inputDirs.length > 1) {
      throw new IllegalArgumentException("Only 1 input directory expected.");
    }

    final File inputDir = inputDirs[0];
    if (inputDir == null || outputDir == null) {
      throw new NullPointerException();
    }

    final HtmlReportHelper hrh = new HtmlReportHelper(outputDir, "index");
    hrh.copyResources(ReportUtils.resourceArray("rtg_logo.png", "rtg.css", "blank.png")); //copy resources up here to create the report files sub dir as well

    final StringBuilder menu = new StringBuilder();

    if (generateTree(inputDir, hrh.getResourcesDir())) {
      menu.append("<a href=\"").append(hrh.getResourcesDirName()).append("/tree.html\">Tree</a><br />");
    }
    if (generateSvd(inputDir, hrh)) {
      menu.append("<a href=\"").append(hrh.getResourcesDirName()).append("/matrix.html\">Matrix</a><br /><a href=\"").append(hrh.getResourcesDirName()).append("/svd.html\">SVD</a><br />");
    }

    final HashMap<String, String> replacements = new HashMap<>();
    replacements.put("__TITLE__", "Similarity");
    replacements.put("__BODY_TEXT__", menu.toString());
    replacements.put("__RESOURCE_DIR__", hrh.getResourcesDirName() + "/");
    ReportUtils.writeHtml(ReportUtils.TEMPLATE_DIR + "/default.html", hrh.getReportFile(), replacements);
  }

  @Override
  public void generateReport(ReportType type, File inputDir, File outputDir) throws IOException {
    generateReport(type, new File[] {inputDir}, outputDir);
  }

  private boolean generateTree(File inputDir, File outputDir) throws IOException {
    final File treeFile = new File(inputDir, "closest.xml");
    if (!treeFile.exists()) {
      Diagnostic.warning("Could not find a closest.xml file in " + inputDir);
      return false;
    } else {
      final StringBuilder tree = new StringBuilder();
      try (BufferedReader br = new BufferedReader(new FileReader(treeFile))) {
        String line;
        while ((line = br.readLine()) != null) {
          line = line.replace("\"", "\\\"");
          tree.append(line);
        }
      }

      // write html
      final HashMap<String, String> replacements = new HashMap<>();
      replacements.put("__TITLE__", "Tree");
      replacements.put("__DATA__", tree.toString());
      ReportUtils.writeHtml(ReportUtils.TEMPLATE_DIR + "/wol.html", new File(outputDir, "tree.html"), replacements);
    }
    return true;
  }

  private boolean generateSvd(File inputDir, HtmlReportHelper hrh) throws IOException {
    final File tsvFile = new File(inputDir, "similarity.tsv");
    if (!tsvFile.exists()) {
      Diagnostic.warning("Could not find a similarity.tsv file in " + inputDir);
      return false;
    } else {
      int numSamples = -1;
      final ArrayList<Heats> heats = new ArrayList<>();
      SimilaritySvd simMatrix = null;

      int min = Integer.MAX_VALUE;
      int max = 0;
      try (BufferedReader br = new BufferedReader(new FileReader(tsvFile))) {
        String line;
        int row = -1;
        while ((line = br.readLine()) != null) {
          if (!line.startsWith("#")) {
            final String[] parts = line.split("\t");
            if (numSamples == -1) {
              numSamples = parts.length - 1;
              simMatrix = new SimilaritySvd(numSamples);
              row = -1;
            } else {
              final Heats h = new Heats(parts[0], numSamples);
              simMatrix.putName(row, parts[0]);
              heats.add(h);
              for (int i = 1; i < parts.length; i++) {
                final int value = Integer.parseInt(parts[i]);
                if (value < min) {
                  min = value;
                }
                if (value > max) {
                  max = value;
                }
                h.setValue(i - 1, value);
                simMatrix.put(i - 1, row, value);
              }
            }
            row++;
          }
        }
      }

      // make min counts middle of colour range
      min -= max - min;
      if (min < 0) {
        min = 0;
      }

      // calculate SVD
      simMatrix.decompose(3);

      //  create plot
      final double[] mins = new double[3];
      final double[] maxs = new double[3];
      for (int i = 0; i < simMatrix.getSvdRowsLength(); i++) {
        for (int j = 0; j < 3; j++) {
          if (i == 0) {
            maxs[j] = mins[j] = simMatrix.getSvdValue(0, j);
          } else {
            final double v = simMatrix.getSvdValue(i, j);
            if (v < mins[j]) {
              mins[j] = v;
            }
            if (v > maxs[j]) {
              maxs[j] = v;
            }
          }
        }
      }
      for (int j = 0; j < 3; j++) {
        if (Double.doubleToLongBits(maxs[j]) == Double.doubleToLongBits(mins[j])) {
          maxs[j] += 1.0;
          mins[j] -= 1.0;
        }
      }

      final StringBuilder samples = new StringBuilder();
      for (int i = 0; i < simMatrix.getSvdRowsLength(); i++) {
        if (i > 0) {
          samples.append(", ");
        }
        samples.append("\"").append(simMatrix.getSvdName(i)).append("\"");
      }
      final StringBuilder data = new StringBuilder();
      for (int i = 0; i < simMatrix.getSvdRowsLength(); i++) {
        if (i > 0) {
          data.append(", ");
        }
        data.append("[");
        for (int j = 0; j < 3; j++) {
          final double v = (simMatrix.getSvdValue(i, j) - mins[j]) / (maxs[j] - mins[j]);
          data.append(String.format(Locale.ROOT, "%1.4f, ", v));
        }
        data.append("\"").append(simMatrix.getSvdName(i)).append("\"]");
      }

      final StringBuilder table = new StringBuilder();
      table.append("<table>").append(StringUtils.LS).append("<tr><th></th>");
      for (int i = 0; i < numSamples; i++) {
        final String name = heats.get(i).getName();
        final String truncName = ReportUtils.truncate(name, 8);
        table.append("<th><a title=\"").append(name).append("\">").append(truncName).append("</a></th>");
      }
      table.append("</tr>");
      for (Heats h : heats) {
        table.append("<tr>").append("<td>").append(h.getName()).append("</td>");
        for (int i = 0; i < numSamples; i++) {
          final String hex = String.format(Locale.ROOT, "%1$02x", 255 - (255 * (h.getValue(i) - min)) / (max - min > 0 ? max - min : 1));
          final String hex2 = String.format(Locale.ROOT, "%1$02x", 255 * (h.getValue(i) - min) / (max - min > 0 ? max - min : 1));
          table.append("<td style=\"background-color: #").append(hex).append(hex2)
          .append("00\"><img src=\"blank.png\" width=\"100\" height=\"15\" title=\"")
          .append(h.getValue(i)).append("\" /></td>");
        }
        table.append("</tr>").append(StringUtils.LS);
      }
      table.append("</table>");

      // write html
      final HashMap<String, String> replacements = new HashMap<>();
      replacements.put("__TITLE__", "Matrix");
      replacements.put("__BODY_TEXT__", table.toString());
      replacements.put("__FILES_DIR__", ".");
      ReportUtils.writeHtml(ReportUtils.TEMPLATE_DIR + "/default.html", new File(hrh.getResourcesDir(), "matrix.html"), replacements);

      replacements.clear();
      replacements.put("__TITLE__", "SVD");
      replacements.put("__SAMPLES__", samples.toString());
      replacements.put("__DATA__", data.toString());
      ReportUtils.writeHtml(ReportUtils.TEMPLATE_DIR + "/svd.html", new File(hrh.getResourcesDir(), "svd.html"), replacements);
      hrh.copyResources(ReportUtils.TEMPLATE_DIR + "/canvasXpress.min.js", ReportUtils.TEMPLATE_DIR + "/excanvas.js");
    }
    return true;
  }
}
