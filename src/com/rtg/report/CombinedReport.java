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
import java.net.URI;
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
 */
public class CombinedReport {

  final List<File> mDirectories;
  final File mOutput;

  /**
   * Create a report for the provided directories
   * @param directories some number of RTG result directories
   * @param output report destination
   */
  public CombinedReport(List<File> directories, File output) {
    mDirectories = directories;
    mOutput = output;
  }

  private static class Filter implements FilenameFilter {
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
      for (String name : FileUtils.list(directory, new Filter())) {
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
    final MetagenomicsReportTemplate template = new MetagenomicsReportTemplate();
    template.mVersion = version;
    template.mMapf = checklist(mapf.size() > 0);
    template.mMapfReport = mapf.size() > 0 ? mapfReport(mapf) : "";
    template.mMap = checklist(map.size() > 0);
    template.mSpecies = checklist(species.size() > 0);
    template.mSpeciesReport = species.size() > 0 ? speciesWrapper(species, map) : "";
    template.mMapx = checklist(mapx.size() > 0);
    template.mMapxReport = mapx.size() > 0 ? mapxReport(mapx) : "";

    final String populated = template.fillTemplate();
    body.append(populated);


    final HtmlReportHelper helper = new HtmlReportHelper(mOutput, "index");
    helper.copyResources(ReportUtils.resourceArray("rtg_logo.png", "rtg.css", "check.png", "cross.png"));
    // write html
    final HashMap<String, String> replace = new HashMap<>();
    replace.put("__TITLE__", "Metagenomics Data Summary Report");
    replace.put("__BODY_TEXT__", body.toString());
    replace.put("__RESOURCE_DIR__", helper.getResourcesDirName() + "/");

    ReportUtils.writeHtml(ReportUtils.TEMPLATE_DIR + "/default.html", helper.getReportFile(), replace);
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
    if (checked) {
      return "checked";
    } else {
      return "crossed";
    }
  }

  private String mapxReport(List<File> mapxDirs) throws IOException {
    final StringBuilder sb = new StringBuilder();
    try (final InputStream resourceAsStream = Resources.getResourceAsStream(ReportUtils.TEMPLATE_DIR + "/mapx_wrapper.html")) {
      final String header = FileUtils.streamToString(resourceAsStream);
      sb.append(header);
    }

    int resultId = 1;
    for (File mapx: mapxDirs) {
      boolean isGzip = false;
      File resultsFile = new File(mapx, ProteinOutputProcessor.TABULAR_ALIGNMENTS);
      if (!resultsFile.exists()) {
        isGzip = true;
       resultsFile = new File(mapx, ProteinOutputProcessor.TABULAR_ALIGNMENTS + FileUtils.GZ_SUFFIX);
      }

      final String fileNameCounter = mapxDirs.size() > 1 ? String.valueOf(resultId) : "";
      final String fileNameSuffix = fileNameCounter + ".tsv" + (isGzip ? ".gz" : "");
      if (!resultsFile.exists()) {
        sb.append(reportSummary(mapx, "Translated search report", "mapxReport" + fileNameCounter, "Mapping summary", null, null, null));
      } else {
        sb.append(reportSummary(mapx, "Translated search report", "mapxReport" + fileNameCounter, "Mapping summary", "Mapping Results", "mapxResults" + fileNameSuffix, resultsFile));
      }
      sb.append("<br/>");
      ++resultId;
    }
    return sb.toString();
  }
  private String mapfReport(List<File> mapfDirs) throws IOException {
    final StringBuilder sb = new StringBuilder();
    sb.append("<h3> Contaminant Removal </h3>\n");
    for (File mapf: mapfDirs) {
      sb.append(reportSummary(mapf, "Filter Report", "mapfReport", "Filter Summary", null, null, null));
    }
    return sb.toString();

  }
  private String reportSummary(File map, String fullName, String fullDirName, String summaryTitle, String resultsLink, String resultFileName, File resultsFile) throws IOException {
    final File relativePath = fullReport(map, fullDirName);
    final String paramsSummary = extractParams(map);
    final MapReportTemplate template = new MapReportTemplate();
    template.mParameterSummary = paramsSummary;
    template.mSummaryFile = summaryFile(map);
    template.mSummaryFileTitle = summaryTitle;
    template.mFullName = fullName;
    final File file = new File(relativePath, "index.html");
    final URI context = mOutput.toURI();
    final URI absolute = file.toURI();
    final URI relative = context.relativize(absolute);
    template.mFullReport = relative.toString();
    final StringBuilder sb = new StringBuilder();
    sb.append(template.fillTemplate());
    if (resultsFile != null) {
      copyFile(resultsFile, new File(mOutput, resultFileName));
      sb.append("<a href=\"")
          .append(resultFileName)
          .append("\">")
          .append(resultsLink)
          .append("</a>");
    }
    return sb.toString();
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

  private File fullReport(File map, String fullDirName) throws IOException {
    final HtmlReportHelper helper = new HtmlReportHelper(map, MapReport.REPORT_NAME);
    final File mapReport1 = new File(mOutput, fullDirName);
    makeOrThrow(mapReport1);
    if (helper.getReportFile().exists() && helper.getResourcesDir().exists()) {
      final HtmlReportHelper destinationReport = new HtmlReportHelper(mapReport1, MapReport.REPORT_NAME);
      copyFile(helper.getReportFile(), destinationReport.getReportFile());
      copyDirRecursive(helper.getResourcesDir(), destinationReport.getResourcesDir());
    } else {
      final MapReport mapReport = new MapReport();
      mapReport.generateReport(ReportType.HTML, mapReport1, map);
    }
    return mapReport1;
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
  private String speciesWrapper(List<File> speciesDirs, List<File> mapDirs) throws IOException {
    final StringBuilder mapReport = new StringBuilder();
    for (File map : mapDirs) {
      mapReport.append(reportSummary(map, "Map Report", "mapReport", "Mapping Summary", null, null, null));
    }
    final StringBuilder speciesReport = new StringBuilder();
    int i = 1;
    for (File species : speciesDirs) {
      final String speciesFileName = "species" + (speciesDirs.size() > 1 ? i : "") + ".tsv";
      speciesReport.append(reportSummary(species, "Species Report", "speciesReport", "Diversity Metrics", "Species Results", speciesFileName, new File(species, SpeciesReport.SPECIES_TSV)));
      ++i;
    }
    final HashMap<String, String> replacements = new HashMap<>();
    replacements.put("__MAP_REPORT__", mapReport.toString());
    replacements.put("__SPECIES_REPORT__", speciesReport.toString());
    return template(replacements, "species_wrapper.html");
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
