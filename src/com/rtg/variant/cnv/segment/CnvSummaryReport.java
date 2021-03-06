/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.variant.cnv.segment;

import static com.rtg.variant.cnv.segment.SegmentVcfOutputFormatter.FORMAT_LOGR;
import static com.rtg.vcf.VcfUtils.INFO_END;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import com.rtg.bed.BedRecord;
import com.rtg.bed.BedWriter;
import com.rtg.util.Environment;
import com.rtg.util.MultiMap;
import com.rtg.util.Utils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.intervals.SequenceNameLocus;
import com.rtg.util.intervals.SequenceNameLocusSimple;
import com.rtg.util.io.FileUtils;
import com.rtg.variant.cnv.CnaType;
import com.rtg.variant.cnv.CnvRecordFilter;
import com.rtg.vcf.AssertVcfSorted;
import com.rtg.vcf.PassOnlyFilter;
import com.rtg.vcf.VcfFilter;
import com.rtg.vcf.VcfFilterIterator;
import com.rtg.vcf.VcfIterator;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;

/**
 * Create a summary report of the amplification status of genes of interest with respect a CNV call set.
 */
class CnvSummaryReport {

  static class GeneStatus {
    final CnaType mType;
    final double mLogR;
    final SequenceNameLocus mSpan;    // Tight span
    final SequenceNameLocus mExtent;  // Maximum span allowed by associated confidence interval
    GeneStatus(CnaType type, double logR, SequenceNameLocus span, SequenceNameLocus extent) {
      mType = type;
      mLogR = logR;
      mSpan = span;
      mExtent = extent;
    }

    @Override
    public String toString() {
      return mType.name() + "," + mLogR + "," + mSpan;
    }
  }

  private static final String VERSION_STRING = "#Version " + Environment.getVersion();
  private static final String CNV_SUMMARY_OUTPUT_VERSION = "v1.1";

  private final ReferenceRanges<String> mRegions;
  private final Map<String, SequenceNameLocus> mByName = new TreeMap<>();
  private final boolean mAlterationsOnly;
  private final double mThreshold;

  CnvSummaryReport(ReferenceRanges<String> regions) {
    this(regions, 0, true);
  }

  CnvSummaryReport(ReferenceRanges<String> regions, double threshold, boolean alterationsOnly) {
    mRegions = regions;
    // Fill in sorted list of gene names
    for (String chr : regions.sequenceNames()) {
      final RangeList<String> rr = regions.get(chr);
      for (RangeList.RangeView<String> region : rr.getRangeList()) {
        if (region.getEnclosingRanges().size() != 1) {
          throw new RuntimeException("Overlapping report regions not supported!");
        }
        final String name = region.getEnclosingRanges().get(0).getMeta();
        if (mByName.containsKey(name)) {
          throw new RuntimeException("Duplicate report regions with name: " + name);
        }
        mByName.put(name, new SequenceNameLocusSimple(chr, region.getStart(), region.getEnd()));
      }
    }
    mThreshold = threshold;
    mAlterationsOnly = alterationsOnly;
  }

  private void writeBedHeader(final BedWriter bw) throws IOException {
    bw.writeln(VERSION_STRING + ", CNV summary BED output " + CNV_SUMMARY_OUTPUT_VERSION);
    if (CommandLine.getCommandLine() != null) {
      bw.writeComment("CL\t" + CommandLine.getCommandLine());
    }
    bw.writeComment("RUN-ID\t" + CommandLine.getRunId());
  }

  /**
   * Summarize interactions between a CNV VCF and the gene regions
   * @param vcfFile input CNV VCF
   * @param output destination file for summary report
   * @throws IOException if an I/O problem occurs
   */
  void report(File vcfFile, File output) throws IOException {
    final MultiMap<String, GeneStatus> geneStatus = new MultiMap<>();
    final List<VcfFilter> svFilt = new ArrayList<>();
    svFilt.add(new AssertVcfSorted());
    svFilt.add(new PassOnlyFilter());
    svFilt.add(new CnvRecordFilter(mRegions.sequenceNames(), false));
    try (VcfIterator vr = new VcfFilterIterator(VcfReader.openVcfReader(vcfFile),  svFilt)) {
      while (vr.hasNext()) {
        final VcfRecord rec = vr.next();

        final CnaType status = CnaType.valueOf(rec);
        final String chr = rec.getSequenceName();
        final int start = rec.getStart() + 1;  // Spec says SV start is the base before the SV.
        final Integer end = VcfUtils.getIntegerInfoFieldFromRecord(rec, INFO_END);

        // Determine intersection
        final RangeList<String> rr = mRegions.get(chr);
        final List<RangeList.RangeView<String>> chrGenes = rr.getFullRangeList();
        for (int hit = rr.findFullRangeIndex(start); hit < chrGenes.size() && chrGenes.get(hit).getStart() < end; hit++) {
          final RangeList.RangeView<String> gene = chrGenes.get(hit);
          if (gene.hasRanges()) {
            final double logR = VcfUtils.getDoubleFormatFieldFromRecord(rec, 0, FORMAT_LOGR);
            if (Math.abs(logR) >= mThreshold) {
              final int[] cipos = VcfUtils.getConfidenceInterval(rec, VcfUtils.INFO_CIPOS);
              final int[] ciend = VcfUtils.getConfidenceInterval(rec, VcfUtils.INFO_CIEND);
              geneStatus.put(gene.getEnclosingRanges().get(0).getMeta(), new GeneStatus(status, logR,
                new SequenceNameLocusSimple(chr, start, end),
                new SequenceNameLocusSimple(chr, start + (cipos == null ? 0 : cipos[0]), end + (ciend == null ? 0 : ciend[1]))));
            }
          }
        }
      }
    }

    // Create gene summary lines as BED records
    final ArrayList<BedRecord> recs = summarizeAsBed(geneStatus);

    // Write the BED summary
    try (final BedWriter bw = new BedWriter(FileUtils.createOutputStream(output))) {
      writeBedHeader(bw);
      bw.writeComment("chrom\tstart\tend\tname\textent\tstatus\tRDR\tLR\tSQS\tsegments");
      for (final BedRecord r : recs) {
        bw.write(r);
      }
    }
  }

  private String scoreAsString(final double lr) {
    return Double.isNaN(lr) ? "" : Utils.realFormat(lr, 4);
  }

  private int lengthOfOverlap(final SequenceNameLocus gene, final SequenceNameLocus call) {
    final int left = Math.max(gene.getStart(), call.getStart());
    final int right = Math.min(gene.getEnd(), call.getEnd());
    return right - left;
  }

  private boolean isFullExtent(final SequenceNameLocus gene, final Collection<GeneStatus> alterations) {
    // Assumes alterations are sorted by overall chromosome position
    int pos = gene.getStart();
    CnaType type = null;
    for (final GeneStatus alteration : alterations) {
      if (type == null) {
        type = alteration.mType;
      } else if (type != alteration.mType) {
        return false; // Mixed DUP/DEL
      }
      final int start = alteration.mExtent.getStart();
      final int end = alteration.mExtent.getEnd();
      if (start > pos) {
        return false; // There is a gap
      }
      pos = end;
    }
    return pos >= gene.getEnd();
  }

  private ArrayList<BedRecord> summarizeAsBed(MultiMap<String, GeneStatus> geneStatus) {
    final ArrayList<BedRecord> recs = new ArrayList<>();
    for (Map.Entry<String, SequenceNameLocus> gene : mByName.entrySet()) {
      final String name = gene.getKey();
      final Collection<GeneStatus> alterations = geneStatus.get(name);
      if (mAlterationsOnly && alterations == null) {
        continue;
      }
      final double lr;
      final SequenceNameLocus region = gene.getValue();
      final String extent;
      final String status;
      final StringBuilder segments = new StringBuilder();
      if (alterations == null) {
        extent = "Full";
        status = "No Alteration";
        lr = Double.NaN;
      } else {
        final boolean gain = alterations.stream().anyMatch(cna -> cna.mType == CnaType.DUP);
        final boolean loss = alterations.stream().anyMatch(cna -> cna.mType == CnaType.DEL);
        extent = isFullExtent(region, alterations) ? "Full" : "Partial";
        status = gain ? loss ? "Mixed" : "Amplification" : "Deletion";
        double best = Double.POSITIVE_INFINITY;
        long bestOverlap = 0;
        for (final GeneStatus alteration : alterations) {
          // Report each call overlapping the gene
          segments.append(alteration).append("; ");
          // Track the best overlap for reporting a final gene-level score
          final long overlap = lengthOfOverlap(gene.getValue(), alteration.mSpan);
          if (overlap >= bestOverlap) {
            if (overlap == bestOverlap) {
              // If tied for length, choose the least modified call
              final double absBest = Math.abs(best);
              final double absCall = Math.abs(alteration.mLogR);
              if (absCall < absBest) {
                best = alteration.mLogR;
              }
            } else {
              bestOverlap = overlap;
              best = alteration.mLogR;
            }
          }
        }
        lr = best;
      }
      recs.add(new BedRecord(region.getSequenceName(), region.getStart(), region.getEnd(), name, extent, status,
        scoreAsString(SegmentVcfOutputFormatter.rdr(lr)), scoreAsString(lr), scoreAsString(SegmentVcfOutputFormatter.sqs(lr)), segments.toString()));
    }
    return recs;
  }
}
