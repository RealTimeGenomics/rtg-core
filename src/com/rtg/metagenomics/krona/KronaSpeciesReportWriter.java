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
package com.rtg.metagenomics.krona;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Map;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.reeltwo.jumble.annotations.JumbleIgnore;
import com.rtg.metagenomics.SpeciesParams;
import com.rtg.report.MapSummaryReport;
import com.rtg.taxonomy.TaxonNode;
import com.rtg.taxonomy.Taxonomy;
import com.rtg.util.HtmlReportHelper;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.io.FileUtils;

/**
 */
@JumbleIgnore
public class KronaSpeciesReportWriter {

  private final HtmlReportHelper mReportHelper;
  private SpeciesParams mParams = null;

  public void setParams(SpeciesParams params) {
    mParams = params;
  }

  /**
   * Create a report writer for Krona output from the Species module
   * @param reportHelper an <code>HtmlReportHelper</code> containing information on the output location etc
   */
  public KronaSpeciesReportWriter(HtmlReportHelper reportHelper) {
    mReportHelper = reportHelper;
  }

  /**
   * Write a <code>Krona</code> report for the given taxonomy, depth first
   * @param rootTaxonNode the root taxonomy node
   * @param taxonIdToNodeMap data map keyed on taxonomy id
   * @throws IOException if an exception occurs while writing
   */
  public void writeReport(TaxonNode rootTaxonNode, Map<Integer, KronaSpeciesNode> taxonIdToNodeMap) throws IOException {
    try (OutputStream os = FileUtils.createOutputStream(mReportHelper.getReportFile())) {
      try {
        final XMLStreamWriter xmlWriter = XMLOutputFactory.newFactory().createXMLStreamWriter(os);
        try {
          header(xmlWriter, true);

          writeNode(xmlWriter, 0, rootTaxonNode, taxonIdToNodeMap);

          xmlWriter.writeCharacters(StringUtils.LS);
          xmlWriter.writeEndElement(); //krona

          htmlFooter(xmlWriter);

          xmlWriter.writeEndDocument();

          copyResources();
        } finally {
          xmlWriter.close();
        }
      } catch (final XMLStreamException e) {
        throw new IOException(e);
      }
    }
  }

  private void writeNode(XMLStreamWriter xmlWriter, int indentLevel, TaxonNode tn, Map<Integer, KronaSpeciesNode> taxonIdToNodeMap) throws XMLStreamException {
    if (tn == null) {
      return;
    }
    final KronaSpeciesNode node = taxonIdToNodeMap.get(tn.getId());
    if (node == null) {
      return;
    }

    xmlWriter.writeCharacters(StringUtils.LS);
    for (int i = 0; i < indentLevel; ++i) {
      xmlWriter.writeCharacters("  ");
    }
    xmlWriter.writeStartElement("node");
    if (tn.getId() == Taxonomy.ROOT_ID) {
      xmlWriter.writeAttribute("name", "all");
    } else {
      xmlWriter.writeAttribute("name", tn.getName());
    }

    //magnitude
    xmlWriter.writeStartElement("magnitude");
    if (node.mAbundance != null) {
      writeXmlVal(xmlWriter, node.mAbundance.toString(), null);
    }
    if (node.mDnaFraction != null) {
      writeXmlVal(xmlWriter, node.mDnaFraction.toString(), null);
    }
    xmlWriter.writeEndElement(); //magnitude

    if (tn.getId() != Taxonomy.ROOT_ID) {
      //taxon id + ncbi link
      xmlWriter.writeStartElement("taxonomyid");
      writeXmlVal(xmlWriter, String.valueOf(tn.getId()), String.valueOf(tn.getId()));
      xmlWriter.writeEndElement(); //taxonomy id

      // rank
      xmlWriter.writeStartElement("rank");
      writeXmlVal(xmlWriter, tn.getRank(), null);
      xmlWriter.writeEndElement(); //abundance

      //abundance
      if (node.mAbundance != null) {
        xmlWriter.writeStartElement("abundance");
        final String abundanceStr = formatDouble(node.mAbundance.getValue())
                  + (tn.getId() == Taxonomy.ROOT_ID ? "" : " (" + formatDouble(node.mAbundance.getLow()) + '~' + formatDouble(node.mAbundance.getHigh()) + ')');
        writeXmlVal(xmlWriter, abundanceStr, null);
        xmlWriter.writeEndElement(); //abundance
      }

      //confidence
      if (node.mConfidence != null) {
        xmlWriter.writeStartElement("confidence");
        writeXmlVal(xmlWriter, Utils.realFormat(node.mConfidence, 2), null);
        xmlWriter.writeEndElement(); //confidence
      }

      //dna fraction
      if (node.mDnaFraction != null) {
        xmlWriter.writeStartElement("dnafraction");
        final String dnafraction = formatDouble(node.mDnaFraction.getValue())
                  + " (" + formatDouble(node.mDnaFraction.getLow()) + '~' + formatDouble(node.mDnaFraction.getHigh()) + ')';
        writeXmlVal(xmlWriter, dnafraction, null);
        xmlWriter.writeEndElement(); //dna fraction
      }

      //coverage depth
      if (node.mCoverageDepth != null) {
        xmlWriter.writeStartElement("coveragedepth");
        writeXmlVal(xmlWriter, formatDouble(node.mCoverageDepth), null);
        xmlWriter.writeEndElement(); //cov depth
      }

      //coverage breadth
      if (node.mCoverageBreadth != null) {
        xmlWriter.writeStartElement("coveragebreadth");
        writeXmlVal(xmlWriter, formatDouble(node.mCoverageBreadth), null);
        xmlWriter.writeEndElement(); //cov breadth
      }

      //genome length
      if (node.mGenomeLength != null) {
        xmlWriter.writeStartElement("genomelength");
        writeXmlVal(xmlWriter, formatLong(node.mGenomeLength), null);
        xmlWriter.writeEndElement(); //genome length
      }

      //mapped reads
      if (node.mMappedReads != null) {
        xmlWriter.writeStartElement("mappedreads");
        writeXmlVal(xmlWriter, formatDouble(node.mMappedReads), null);
        xmlWriter.writeEndElement(); //mapped reads
      }
    }
    for (final TaxonNode child : tn.getChildren()) {
      writeNode(xmlWriter, indentLevel + 1, child, taxonIdToNodeMap);
    }

    xmlWriter.writeEndElement(); //node
  }

  private String formatDouble(Double d) {
    if (d == null) {  //this shouldn't ever happen...
      return "0.00";
    }
    return Utils.realFormat(d, 2);
  }

  private String formatLong(Long d) {
    if (d == null) {  //this shouldn't ever happen...
      return "0";
    }
    return d.toString();
  }

  void writeXmlVal(XMLStreamWriter xmlWriter, String val, String href) throws XMLStreamException {
    xmlWriter.writeStartElement("val");
    if (href != null) {
      xmlWriter.writeAttribute("href", href);
    }
    xmlWriter.writeCharacters(val);
    xmlWriter.writeEndElement();
  }


  private void header(XMLStreamWriter xmlWriter, boolean htmlOutput) throws XMLStreamException {
    if (htmlOutput) {
      htmlHeader(xmlWriter);
    } else {
      xmlWriter.writeStartDocument();
    }
    xmlWriter.writeStartElement("krona");
    xmlWriter.writeCharacters(StringUtils.LS);
    xmlWriter.writeStartElement("attributes");

    xmlWriter.writeAttribute("magnitude", "magnitude");

    xmlWriter.writeStartElement("list");
    xmlWriter.writeCharacters("members");
    xmlWriter.writeEndElement(); //list

    xmlWriter.writeCharacters(StringUtils.LS);
    xmlWriter.writeStartElement("attribute");
    xmlWriter.writeCharacters("magnitude");
    xmlWriter.writeEndElement(); //magnitude attr

    xmlWriter.writeCharacters(StringUtils.LS);
    xmlWriter.writeStartElement("attribute");
    xmlWriter.writeAttribute("display", "Taxonomy Id");
    xmlWriter.writeAttribute("hrefBase", "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=");
    xmlWriter.writeAttribute("mono", "true");
    xmlWriter.writeCharacters("taxonomyid");
    xmlWriter.writeEndElement(); //taxon id attr

    xmlWriter.writeCharacters(StringUtils.LS);
    xmlWriter.writeStartElement("attribute");
    xmlWriter.writeAttribute("display", "Rank");
    xmlWriter.writeAttribute("mono", "true");
    xmlWriter.writeCharacters("rank");
    xmlWriter.writeEndElement(); //rank attr

    xmlWriter.writeCharacters(StringUtils.LS);
    xmlWriter.writeStartElement("attribute");
    xmlWriter.writeAttribute("display", "Abundance");
    xmlWriter.writeAttribute("mono", "true");
    xmlWriter.writeCharacters("abundance");
    xmlWriter.writeEndElement(); //abundance attr

    xmlWriter.writeCharacters(StringUtils.LS);
    xmlWriter.writeStartElement("attribute");
    xmlWriter.writeAttribute("display", "Confidence");
    xmlWriter.writeAttribute("mono", "true");
    xmlWriter.writeCharacters("confidence");
    xmlWriter.writeEndElement(); //confidence

    xmlWriter.writeCharacters(StringUtils.LS);
    xmlWriter.writeStartElement("attribute");
    xmlWriter.writeAttribute("display", "DNA Fraction");
    xmlWriter.writeAttribute("mono", "true");
    xmlWriter.writeCharacters("dnafraction");
    xmlWriter.writeEndElement(); //dnafraction attr

    xmlWriter.writeCharacters(StringUtils.LS);
    xmlWriter.writeStartElement("attribute");
    xmlWriter.writeAttribute("display", "Coverage Depth");
    xmlWriter.writeAttribute("mono", "true");
    xmlWriter.writeCharacters("coveragedepth");
    xmlWriter.writeEndElement(); //coverage depth

    xmlWriter.writeCharacters(StringUtils.LS);
    xmlWriter.writeStartElement("attribute");
    xmlWriter.writeAttribute("display", "Coverage Breadth");
    xmlWriter.writeAttribute("mono", "true");
    xmlWriter.writeCharacters("coveragebreadth");
    xmlWriter.writeEndElement(); //coverage breadth

    xmlWriter.writeCharacters(StringUtils.LS);
    xmlWriter.writeStartElement("attribute");
    xmlWriter.writeAttribute("display", "Genome Length");
    xmlWriter.writeAttribute("mono", "true");
    xmlWriter.writeCharacters("genomelength");
    xmlWriter.writeEndElement(); //genome length

    xmlWriter.writeCharacters(StringUtils.LS);
    xmlWriter.writeStartElement("attribute");
    xmlWriter.writeAttribute("display", "Mapped Reads");
    xmlWriter.writeAttribute("mono", "true");
    xmlWriter.writeCharacters("mappedreads");
    xmlWriter.writeEndElement(); //mappedreads

    xmlWriter.writeCharacters(StringUtils.LS);
    xmlWriter.writeEndElement(); //attrs

    xmlWriter.writeCharacters(StringUtils.LS);
    xmlWriter.writeStartElement("datasets");
    xmlWriter.writeStartElement("dataset");
    xmlWriter.writeCharacters("by abundance");
    xmlWriter.writeEndElement(); //dataset;

    xmlWriter.writeStartElement("dataset");
    xmlWriter.writeCharacters("by DNA fraction");
    xmlWriter.writeEndElement(); //dataset

    xmlWriter.writeEndElement(); //datasets
  }

  private void htmlHeader(XMLStreamWriter writer) throws XMLStreamException {
    writer.writeDTD("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">");
    writer.writeCharacters(StringUtils.LS);

    writer.writeStartElement("html");
    writer.writeAttribute("xmlns", "http://www.w3.org/1999/xhtml");
    writer.writeAttribute("xml:lang", "en");
    writer.writeAttribute("lang", "en");
    writer.writeCharacters(StringUtils.LS);

    writer.writeStartElement("head");
    writer.writeCharacters(StringUtils.LS);

    writer.writeEmptyElement("meta");
    writer.writeAttribute("charset", "utf-8");
    writer.writeCharacters(StringUtils.LS);

    writer.writeEmptyElement("link");
    writer.writeAttribute("rel", "shortcut icon");
    writer.writeAttribute("href", "http://realtimegenomics.com/asset/icon.ico");
    writer.writeCharacters(StringUtils.LS);

    writer.writeEmptyElement("link");
    writer.writeAttribute("rel", "stylesheet");
    writer.writeAttribute("type", "text/css");
    writer.writeAttribute("href", mReportHelper.getResourcesDirName() + "/rtg.css");
    writer.writeCharacters(StringUtils.LS);

    writer.writeStartElement("script");
    writer.writeAttribute("id", "notfound");
    writer.writeCharacters("window.onload=function(){document.body.innerHTML=\"Could not get resource krona-2.0.js\"}");
    writer.writeEndElement(); //script
    writer.writeCharacters(StringUtils.LS);

    writer.writeStartElement("script");
    writer.writeAttribute("src", mReportHelper.getResourcesDirName() + "/krona-2.0.js");
    writer.writeEndElement(); //script
    writer.writeCharacters(StringUtils.LS);

    writer.writeEndElement(); //head
    writer.writeCharacters(StringUtils.LS);
    writer.writeStartElement("body");
    writer.writeCharacters(StringUtils.LS);

    writer.writeStartElement("div");
    writer.writeAttribute("id", "dheader");
    writer.writeStartElement("a");
    writer.writeAttribute("id", "lheader");
    writer.writeAttribute("href", "http://www.realtimegenomics.com");
    writer.writeStartElement("h1");
    writer.writeCharacters("Species Report");
    writer.writeEndElement(); //h1
    writer.writeEndElement(); //a
    writer.writeEndElement(); //div

    writeParams(writer);

    writer.writeEmptyElement("img");
    writer.writeAttribute("id", "hiddenImage");
    writer.writeAttribute("src", mReportHelper.getResourcesDirName() + "/hidden.png");
    writer.writeAttribute("style", "display:none");
    writer.writeCharacters(StringUtils.LS);

    writer.writeEmptyElement("img");
    writer.writeAttribute("id", "loadingImage");
    writer.writeAttribute("src", mReportHelper.getResourcesDirName() + "/loading.gif");
    writer.writeAttribute("style", "display:none");
    writer.writeCharacters(StringUtils.LS);

    writer.writeStartElement("noscript");
    writer.writeCharacters("Javascript must be enabled to view this page.");
    writer.writeEndElement();
    writer.writeCharacters(StringUtils.LS);

    writer.writeStartElement("div");
    writer.writeAttribute("style", "display:none");
    writer.writeCharacters(StringUtils.LS);
  }

  private void writeParams(XMLStreamWriter writer) throws XMLStreamException {
    if (mParams != null) {
      writer.writeStartElement("div");
      writer.writeAttribute("style", "margin:1em;");
      // tag used by combined report to extract parameter summary
      writer.writeCharacters(StringUtils.LS);
      writer.writeComment(MapSummaryReport.SUMMARY_STRING);
      writer.writeCharacters(StringUtils.LS);

      writer.writeStartElement("strong");
      writer.writeCharacters("Sample: ");
      writer.writeEndElement(); //strong
      String join = "";
      for (final File f : mParams.mapped()) {
        writer.writeCharacters(join);
        writer.writeCharacters(f.getPath());
        join = ", ";
      }
      writer.writeEmptyElement("br");
      writer.writeStartElement("strong");
      writer.writeCharacters("Results directory: ");
      writer.writeEndElement(); //strong
      writer.writeCharacters(mParams.directory().getPath());
      writer.writeEmptyElement("br");

      writer.writeStartElement("p");
      writer.writeStartElement("strong");
      writer.writeCharacters("Command line: ");
      writer.writeEndElement(); //strong
      writer.writeCharacters(CommandLine.getCommandLine());
      writer.writeEndElement(); //p

      // tag used by combined report to extract parameter summary
      writer.writeCharacters(StringUtils.LS);
      writer.writeComment(MapSummaryReport.END_SUMMARY_STRING);
      writer.writeCharacters(StringUtils.LS);
      writer.writeEndElement(); // div
    }
  }

  private void htmlFooter(XMLStreamWriter writer) throws XMLStreamException {
    writer.writeEndElement(); //div



    writer.writeEndElement(); //body
    writer.writeEndElement(); //html
  }

  void copyResources() throws IOException {
    mReportHelper.copyResources("com/rtg/report/resources/rtg.css",
                                      "com/rtg/report/resources/rtg_logo.png",
                                      "com/rtg/metagenomics/krona/resources/krona-2.0.js",
                                      "com/rtg/metagenomics/krona/resources/hidden.png",
                                      "com/rtg/metagenomics/krona/resources/loading.gif");
  }
}
