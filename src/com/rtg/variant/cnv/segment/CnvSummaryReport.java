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
import static com.rtg.variant.cnv.segment.SegmentVcfOutputFormatter.INFO_END;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import com.rtg.bed.BedRecord;
import com.rtg.bed.BedWriter;
import com.rtg.sam.SamRangeUtils;
import com.rtg.util.MultiMap;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.intervals.SequenceNameLocus;
import com.rtg.util.intervals.SequenceNameLocusSimple;
import com.rtg.util.io.FileUtils;
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
    final boolean mPartial;
    final double mLogR;
    final SequenceNameLocus mSpan;
    GeneStatus(CnaType type, double logR, boolean partial, SequenceNameLocus span) {
      mType = type;
      mLogR = logR;
      mPartial = partial;
      mSpan = span;
    }
    public String toString() {
      return (mPartial ? "Partial" : "Full") + "," + mType.name() + "," + mLogR + "," + mSpan;
    }
  }


  private final ReferenceRanges<String> mRegions;
  private final Map<String, SequenceNameLocus> mByName = new TreeMap<>();
  private final double mThreshold;

  CnvSummaryReport(ReferenceRanges<String> regions) {
    this(regions, 0);
  }

  CnvSummaryReport(ReferenceRanges<String> regions, double threshold) {
    mRegions = regions;
    // Fill in sorted list of gene names
    for (String chr : regions.sequenceNames()) {
      final RangeList<String> rr = regions.get(chr);
      for (RangeList.RangeData<String> region : rr.getRangeList()) {
        if (region.getMeta().size() != 1) {
          throw new RuntimeException("Overlapping report regions not supported!");
        }
        final String name = region.getMeta().get(0);
        if (mByName.containsKey(name)) {
          throw new RuntimeException("Duplicate report regions with name: " + name);
        }
        mByName.put(name, new SequenceNameLocusSimple(chr, region.getStart(), region.getEnd()));
      }
    }
    mThreshold = threshold;
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
    svFilt.add(new CnvRecordFilter(mRegions.sequenceNames()));
    try (VcfIterator vr = new VcfFilterIterator(VcfReader.openVcfReader(vcfFile),  svFilt)) {
      while (vr.hasNext()) {
        final VcfRecord rec = vr.next();

        final CnaType status = CnaType.valueOf(rec);
        final String chr = rec.getSequenceName();
        final int start = rec.getStart() + 1;  // Spec says SV start is the base before the SV.
        final Integer end = VcfUtils.getIntegerInfoFieldFromRecord(rec, INFO_END) - 1; // Convert from 1-based to 0-based

        // Determine intersection
        final RangeList<String> rr = mRegions.get(chr);
        final List<RangeList.RangeData<String>> chrGenes = rr.getFullRangeList();
        for (int hit = rr.findFullRangeIndex(start); hit < chrGenes.size() && chrGenes.get(hit).getStart() < end; hit++) {
          final RangeList.RangeData<String> gene = chrGenes.get(hit);
          if (gene.getMeta() != null) {
            final double logR = VcfUtils.getDoubleFormatFieldFromRecord(rec, 0, FORMAT_LOGR);
            if (Math.abs(logR) >= mThreshold) {
              final boolean partial = gene.getStart() < start || gene.getEnd() > end;
              geneStatus.put(gene.getMeta().get(0), new GeneStatus(status, logR, partial, new SequenceNameLocusSimple(chr, start, end)));
            }
          }
        }
      }
    }

    // Create gene summary lines as BED records
    final ArrayList<BedRecord> recs = summarizeAsBed(geneStatus);

    //recs.sort(SequenceNameLocusComparator.SINGLETON);

    // Write the BED summary
    try (BedWriter bw = new BedWriter(FileUtils.createOutputStream(output))) {
      bw.writeComment("crom\tstart\tend\tgene\tstatus\tsegments");
      for (BedRecord r : recs) {
        bw.write(r);
      }
    }
  }

  private ArrayList<BedRecord> summarizeAsBed(MultiMap<String, GeneStatus> geneStatus) {
    final ArrayList<BedRecord> recs = new ArrayList<>();
    for (Map.Entry<String, SequenceNameLocus> gene : mByName.entrySet()) {
      final String name = gene.getKey();
      final Collection<GeneStatus> alterations = geneStatus.get(name);
      final SequenceNameLocus region = gene.getValue();
      final String status;
      final StringBuilder segments = new StringBuilder();
      if (alterations == null) {
        status = "No Alteration";
      } else {
        final boolean gain = alterations.stream().anyMatch(cna -> cna.mType == CnaType.DUP);
        final boolean loss = alterations.stream().anyMatch(cna -> cna.mType == CnaType.DEL);
        status = gain ? loss ? "Mixed" : "Amplification" : "Deletion";
        for (GeneStatus alteration : alterations) {
          segments.append(alteration).append("; ");
        }
      }
      recs.add(new BedRecord(region.getSequenceName(), region.getStart(), region.getEnd(), name, status, segments.toString()));
    }
    return recs;
  }


  /**
   * Direct entry point
   * @param args arguments
   * @throws IOException if the wheels fall off
   */
  public static void main(final String[] args) throws IOException {
    final File geneBed = new File(args[0]);
    final File vcfFile = new File(args[1]);
    final double threshold = Double.parseDouble(args[2]);
    final CnvSummaryReport p = new CnvSummaryReport(SamRangeUtils.createBedReferenceRanges(geneBed), threshold);
    p.report(vcfFile, new File("-"));
  }
}
