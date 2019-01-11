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
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.rtg.protein.ProteinOutputProcessor;
import com.rtg.util.Environment;
import com.rtg.util.HtmlReportHelper;
import com.rtg.util.Resources;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LineWriter;

/**
 * Create a combined report from several pre-existing output directories.
 * This is extremely brittle.
 */
public class CombinedReport {

  private final List<File> mDirectories;
  private final File mOutput;
  private final HtmlReportHelper mMainHelper;

  /**
   * Create a report for the provided directories
   * @param directories some number of RTG result directories
   * @param output report destination
   */
  public CombinedReport(List<File> directories, File output) {
    mDirectories = directories;
    mOutput = output;
    mMainHelper = new HtmlReportHelper(mOutput, "index");
  }

  private static class LogfileFilter implements FilenameFilter {
    @Override
    public boolean accept(File dir, String name) {
      return name.endsWith(".log");
    }
  }

  /**
   * Perform report writing
   * @throws IOException when one of the myriad of IO operations fails
   */
  public void makeReport() throws IOException {
    final StringBuilder body = new StringBuilder();
    final String version = Environment.getVersion();

    final List<File> mapf = new ArrayList<>();
    final List<File> map = new ArrayList<>();
    final List<File> species = new ArrayList<>();
    final List<File> mapx = new ArrayList<>();
    for (File directory : mDirectories) {
      for (String name : FileUtils.list(directory, new LogfileFilter())) {
        switch(name) {
          case "mapf.log" :
            mapf.add(directory);
            break;
          case "mapx.log" :
            mapx.add(directory);
            break;
          case "map.log" :
            map.add(directory);
            break;
          case "species.log" :
            species.add(directory);
            break;
          default:
            // Can't do anything about directories we don't recognize
            break;
        }
      }
    }
    mMainHelper.copyResources(ReportUtils.resourceArray("rtg.css"));

    final MetagenomicsReportTemplate template = new MetagenomicsReportTemplate();
    template.mVersion = version;
    template.mMapf = checklist(mapf.size() > 0);
    template.mMapfReport = mapf.size() > 0 ? mapfReport(mapf) : "";
    template.mMap = checklist(map.size() > 0);
    template.mSpecies = checklist(species.size() > 0);
    template.mSpeciesReport = species.size() > 0 ? speciesReport(species, map) : "";
    template.mMapx = checklist(mapx.size() > 0);
    template.mMapxReport = mapx.size() > 0 ? mapxReport(mapx) : "";

    final String populated = template.fillTemplate();
    body.append(populated);

    // write html
    final HashMap<String, String> replace = new HashMap<>();
    replace.put("__TITLE__", "Metagenomics Data Summary Report");
    replace.put("__BODY_TEXT__", body.toString());
    replace.put("__RESOURCE_DIR__", mMainHelper.getResourcesDirName() + "/");

    ReportUtils.writeHtml(ReportUtils.TEMPLATE_DIR + "/default.html", mMainHelper.getReportFile(), replace);
  }

  static String template(Map<String, String> replacements, String templateName) throws IOException {
    final StringBuilder body = new StringBuilder();
    final String templateFile = ReportUtils.TEMPLATE_DIR + "/" + templateName;
    try (ByteArrayOutputStream baos = new ByteArrayOutputStream()) {
      try (
        LineWriter writer = new LineWriter(new OutputStreamWriter(baos));
        BufferedReader template = new BufferedReader(new InputStreamReader(Resources.getResourceAsStream(SpeciesReport.class, templateFile)))
      ) {
        ReportUtils.fillTemplate(replacements, template, writer);
      }
      body.append(baos);
    }
    return body.toString();
  }

  private static void makeOrThrow(File dir) throws IOException {
    if (!dir.exists() && !dir.mkdir()) {
      throw new IOException("Can't create directory: '" + dir.getPath() + "'");
    }
  }

  private String checklist(boolean checked) {
    return checked ? "checked" : "crossed";
  }

  static final class ReportResultsFile {
    private final File mSourceFile;
    private final String mDestFilename;
    private final String mDestLinkTitle;
    static ReportResultsFile getResultsFile(File sourceFile, String destFilename, String linkTitle) {
      if (sourceFile.exists()) {
        return new ReportResultsFile(sourceFile, destFilename, linkTitle);
      } else {
        final File newsourceFile = new File(sourceFile.getParent(), sourceFile.getName() + FileUtils.GZ_SUFFIX);
        if (newsourceFile.exists()) {
          return new ReportResultsFile(newsourceFile, destFilename + FileUtils.GZ_SUFFIX, linkTitle);
        }
      }
      return null;
    }
    ReportResultsFile(File sourceFile, String destFilename, String linkTitle) {
        mSourceFile = sourceFile;
        mDestFilename = destFilename;
        mDestLinkTitle = linkTitle;
    }
  }

  private String speciesReport(List<File> speciesDirs, List<File> mapDirs) throws IOException {
    final StringBuilder mapReport = new StringBuilder();
    for (File map : mapDirs) {
      mapReport.append(reportSummary(map, "Map Report", "mapReport", "Mapping Summary"));
    }
    final StringBuilder speciesReport = new StringBuilder();
    int counter = 0;
    for (File species : speciesDirs) {
      final String id = speciesDirs.size() > 1 ? String.valueOf(++counter) : "";
      final ReportResultsFile resultsFile = ReportResultsFile.getResultsFile(new File(species, SpeciesReport.SPECIES_TSV), "species" + id + ".tsv", "Species Results");
      speciesReport.append(reportSummary(species, "Species Report", "speciesReport" + id, "Diversity Metrics", resultsFile));
    }
    final HashMap<String, String> replacements = new HashMap<>();
    replacements.put("__MAP_REPORT__", mapReport.toString());
    replacements.put("__SPECIES_REPORT__", speciesReport.toString());
    return template(replacements, "species_wrapper.html");
  }
  private String mapxReport(List<File> mapxDirs) throws IOException {
    final StringBuilder sb = new StringBuilder();
    try (final InputStream resourceAsStream = Resources.getResourceAsStream(ReportUtils.TEMPLATE_DIR + "/mapx_wrapper.html")) {
      final String header = FileUtils.streamToString(resourceAsStream);
      sb.append(header);
    }
    int counter = 0;
    for (File mapx : mapxDirs) {
      final String id = mapxDirs.size() > 1 ? String.valueOf(++counter) : "";
      final ReportResultsFile resultsFile = ReportResultsFile.getResultsFile(new File(mapx, ProteinOutputProcessor.TABULAR_ALIGNMENTS), "mapxResults" + id + ".tsv", "Mapping Results");
      sb.append(reportSummary(mapx, "Translated search report", "mapxReport" + id, "Mapping summary", resultsFile));
    }
    return sb.toString();
  }
  private String mapfReport(List<File> mapfDirs) throws IOException {
    final StringBuilder sb = new StringBuilder();
    sb.append("<h3> Contaminant Removal </h3>\n");
    int counter = 0;
    for (File mapf: mapfDirs) {
      final String id = mapfDirs.size() > 1 ? String.valueOf(++counter) : "";
      sb.append(reportSummary(mapf, "Filter Report", "mapfReport" + id, "Filter Summary"));
    }
    return sb.toString();
  }

  private String reportSummary(File sourceDir, String reportTitle, String destDirName, String summaryTitle, ReportResultsFile... resultsFiles) throws IOException {
    final File relativePath = fullReport(sourceDir, destDirName);
    final MapReportTemplate template = new MapReportTemplate();
    template.mParameterSummary = extractParams(sourceDir);
    template.mSummaryFileTitle = summaryTitle;
    template.mSummaryFile = summaryFile(sourceDir);
    template.mFullReportTitle = reportTitle;
    template.mFullReport = mOutput.toURI().relativize(new File(relativePath, "index.html").toURI()).toString();
    final StringBuilder sb = new StringBuilder();
    sb.append(template.fillTemplate());
    for (ReportResultsFile resultsFile : resultsFiles) {
      if (resultsFile != null) {
        copyResultsFile(resultsFile, sb);
      }
    }
    return sb.toString();
  }

  private void copyResultsFile(ReportResultsFile resultsFile, StringBuilder htmlOut) throws IOException {
    final File destFile = new File(mMainHelper.getResourcesDir(), resultsFile.mDestFilename);
    copyFile(resultsFile.mSourceFile, destFile);
    htmlOut.append("<a href=\"").append(mOutput.toURI().relativize(destFile.toURI()).toString()).append("\">").append(resultsFile.mDestLinkTitle).append("</a>").append("<br/>");
  }

  private String extractParams(File map) throws IOException {
    final HtmlReportHelper helper = new HtmlReportHelper(map, MapReport.REPORT_NAME);
    if (!helper.getReportFile().exists()) {
      return "";
    } else {
      final StringBuilder sb = new StringBuilder();
      try (BufferedReader reader = new BufferedReader(FileUtils.createReader(helper.getReportFile(), false))) {
        String line;
        while ((line = reader.readLine()) != null && !MapReport.START_PARAMS_SUMMARY.equals(line)) {
          //skip line
        }

        while ((line = reader.readLine()) != null && !MapReport.END_PARAMS_SUMMARY.equals(line)) {
          sb.append(line).append(StringUtils.LS);
        }
      }
      return sb.toString();
    }
  }

  private File fullReport(File sourceDir, String destDirName) throws IOException {
    final File destDir = new File(mMainHelper.getResourcesDir(), destDirName);
    makeOrThrow(destDir);
    final HtmlReportHelper helper = new HtmlReportHelper(sourceDir, MapReport.REPORT_NAME);
    if (helper.getReportFile().exists() && helper.getResourcesDir().exists()) {
      final HtmlReportHelper destinationReport = new HtmlReportHelper(destDir, MapReport.REPORT_NAME);
      copyFile(helper.getReportFile(), destinationReport.getReportFile());
      copyDirRecursive(helper.getResourcesDir(), destinationReport.getResourcesDir());
    }
    return destDir;
  }

  static String getCommandLine(File logFile) throws IOException {
    String commandline = "";
    try (BufferedReader reader = new BufferedReader(new FileReader(logFile))) {
      String line;
      while ((line = reader.readLine()) != null) {
        if (line.contains("Command line arguments:")) {
          commandline = line.replaceAll("^.* Command line arguments: ", "");
          break;
        }
      }
    }
    return commandline;
  }

  String summaryFile(File dir) throws IOException {
    final File summary = new File(dir, "summary.txt");
    if (summary.exists()) {
      return FileUtils.fileToString(summary);
    } else {
      return "";
    }
  }

  private void copyDirRecursive(File src, File dest) throws IOException {
    makeOrThrow(dest);
    final File[] list = src.listFiles();
    if (list != null) {
      for (File f : list) {
        final File destFile = new File(dest, f.getName());
        if (f.isDirectory()) {
          copyDirRecursive(f, destFile);
        } else {
          copyFile(f, destFile);
        }
      }
    }
  }

  private void copyFile(File src, File dest) throws IOException {
    try (InputStream in = new FileInputStream(src);
         OutputStream out = new FileOutputStream(dest)) {
      FileUtils.streamToStream(in, out, 1000);
    }
  }

  /**
   * Command line entry point
   * @param args command line arguments
   * @throws IOException at will
   */
  public static void  main(String[] args) throws IOException {
    final File output = new File(args[0]);
    if (output.exists()) {
      throw new NoTalkbackSlimException("output dir exists");
    } else {
      makeOrThrow(output);
    }
    final List<File> reports = new ArrayList<>(Math.max(0, args.length - 1));
    for (int i = 1; i < args.length; ++i) {
      reports.add(new File(args[i]));
    }
    final CombinedReport combined = new CombinedReport(reports, output);
    combined.makeReport();

  }
}
