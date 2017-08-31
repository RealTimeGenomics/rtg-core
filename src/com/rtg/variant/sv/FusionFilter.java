/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.sv;

import static com.rtg.launcher.CommonFlags.FILE;
import static com.rtg.launcher.CommonFlags.INPUT_FLAG;
import static com.rtg.launcher.CommonFlags.INT;
import static com.rtg.launcher.CommonFlags.NO_GZIP;
import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;
import static com.rtg.vcf.VcfUtils.INFO_SVTYPE;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.rtg.bed.BedReader;
import com.rtg.bed.BedRecord;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.util.Pair;
import com.rtg.util.TsvParser;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.io.FileUtils;
import com.rtg.vcf.AltVariantTypeFilter;
import com.rtg.vcf.AssertVcfSorted;
import com.rtg.vcf.BreakpointAlt;
import com.rtg.vcf.VariantType;
import com.rtg.vcf.VcfAnnotator;
import com.rtg.vcf.VcfFilter;
import com.rtg.vcf.VcfFilterIterator;
import com.rtg.vcf.VcfIterator;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.VcfWriter;
import com.rtg.vcf.VcfWriterFactory;
import com.rtg.vcf.header.InfoField;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

/**
 * Filters structural variant VCF for putative gene fusions based on location and orientation.
 */
public class FusionFilter extends AbstractCli {

  private static final String INFO_ACCEPTOR_GENE = "ACCEPTOR_GENE";
  private static final String INFO_DONOR_GENE = "DONOR_GENE";
  private static final String INFO_FUSION_SCORE = "FUSION_SCORE";
  private static final String INFO_FUSION_NAME = "FUSION_NAME";

  private static final String ACCEPTOR_REGIONS = "acceptor-regions";
  private static final String DONOR_REGIONS = "donor-regions";
  private static final String SCORES = "fusion-scores";

  private static final String NO_HEADER = "no-header";
  private static final String GENE_COLUMN = "gene-column";
  private static final String STRAND_COLUMN = "strand-column";

  @Override
  public String moduleName() {
    return "fusionfilter";
  }

  @Override
  public String description() {
    return "filter a structural variant VCF for putative gene fusions";
  }

  @Override
  protected void initFlags() {
    CommonFlagCategories.setCategories(mFlags);
    mFlags.setDescription("Filters a VCF containing structural variant records to output only those that could correspond to a gene fusion between a donor and acceptor gene.");
    mFlags.registerRequired('i', INPUT_FLAG, File.class, FILE, "VCF file containing variants to filter. Use '-' to read from standard input").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('o', OUTPUT_FLAG, File.class, FILE, "output VCF file name. Use '-' to write to standard output").setCategory(INPUT_OUTPUT);
    CommonFlags.initForce(mFlags);
    mFlags.registerRequired('a', ACCEPTOR_REGIONS, File.class, FILE, "BED file specifying strand-specific acceptor gene regions").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('d', DONOR_REGIONS, File.class, FILE, "BED file specifying strand-specific donor gene regions").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional('s', SCORES, File.class, FILE, "if set, annotate fusions with fusion-specific scores supplied in the supplied file").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(GENE_COLUMN, Integer.class, INT, "column within the BED file that supplies gene name", 4).setCategory(REPORTING);
    mFlags.registerOptional(STRAND_COLUMN, Integer.class, INT, "column within the BED file that supplies strand", 6).setCategory(REPORTING);
    mFlags.registerOptional(NO_HEADER, "prevent VCF header from being written").setCategory(UTILITY);
    CommonFlags.initNoGzip(mFlags);
    mFlags.setValidator(flags -> CommonFlags.validateInputFile(flags, INPUT_FLAG, ACCEPTOR_REGIONS, DONOR_REGIONS, SCORES)
      && CommonFlags.validateOutputFile(flags, VcfUtils.getZippedVcfFileName(!flags.isSet(NO_GZIP), (File) flags.getValue(OUTPUT_FLAG)))
      && mFlags.checkInRange(GENE_COLUMN, 4, Integer.MAX_VALUE)
      && mFlags.checkInRange(STRAND_COLUMN, 4, Integer.MAX_VALUE));
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final File inputFile = (File) mFlags.getValue(INPUT_FLAG);
    final File output = (File) mFlags.getValue(OUTPUT_FLAG);
    final boolean gzip = !mFlags.isSet(NO_GZIP);
    final boolean stdout = FileUtils.isStdio(output);

    final List<VcfFilter> filters = new ArrayList<>();
    filters.add(new AssertVcfSorted());
    filters.add(new AltVariantTypeFilter(EnumSet.of(VariantType.SV_BREAKEND, VariantType.SV_SYMBOLIC)));
    filters.add(new RegionFilter((File) mFlags.getValue(ACCEPTOR_REGIONS), (File) mFlags.getValue(DONOR_REGIONS),
      (Integer) mFlags.getValue(GENE_COLUMN) - 4, (Integer) mFlags.getValue(STRAND_COLUMN) - 4));

    final List<VcfAnnotator> annotators = new ArrayList<>();
    if (mFlags.isSet(SCORES)) {
      annotators.add(new FusionAnnotator((File) mFlags.getValue(SCORES)));
    }

    try (VcfIterator reader = new VcfFilterIterator(VcfReader.openVcfReader(inputFile), filters)) {
      final VcfHeader header = reader.getHeader();
      header.addRunInfo();
      filters.forEach(filter -> filter.setHeader(header));
      annotators.forEach(ann -> ann.updateHeader(header));
      final File vcfFile = stdout ? null : VcfUtils.getZippedVcfFileName(gzip, output);
      try (VcfWriter writer = new VcfWriterFactory(mFlags).addRunInfo(true).make(header, vcfFile, out)) {
        while (reader.hasNext()) {
          final VcfRecord rec = reader.next();
          annotators.forEach(filter -> filter.annotate(rec));
          writer.write(rec);
        }
      }
    }
    return 0;
  }

  // Annotate recognized gene fusions with their name and score
  private static class FusionAnnotator implements VcfAnnotator {
    private static final class ScoreParser extends TsvParser<Map<String, Integer>> {
      private final Map<String, Integer> mMap = new HashMap<>();
      @Override
      protected void parseLine(String... parts) throws IOException {
        if (parts.length != 2) {
          throw new IOException("Expected 2 columns on line: " + line());
        }
        final int score;
        try {
          score = Integer.parseInt(parts[1]);
        } catch (NumberFormatException e) {
          throw new IOException("Malformed score on line: " + line(), e);
        }
        if (mMap.put(parts[0], score) != null) {
          throw new IOException("Duplicate fusion name detected: " + line());
        }
      }
      @Override
      protected Map<String, Integer> result() {
        return mMap;
      }
    }

    private final Map<String, Integer> mFusionScores;
    FusionAnnotator(File scoreFile) throws IOException {
      mFusionScores = new ScoreParser().parse(scoreFile);
    }
    @Override
    public void updateHeader(VcfHeader header) {
      header.ensureContains(new InfoField(INFO_FUSION_SCORE, MetaType.INTEGER, VcfNumber.ONE, "Score assigned to the gene fusion"));
      header.ensureContains(new InfoField(INFO_FUSION_NAME, MetaType.STRING, VcfNumber.ONE, "Name of the highest scoring gene fusion"));
    }

    @Override
    public void annotate(VcfRecord rec) {
      Integer bestScore = 0;
      String bestFusion = null;
      for (String a : rec.getInfo().get(INFO_ACCEPTOR_GENE)) {
        for (String d : rec.getInfo().get(INFO_DONOR_GENE)) {
          final String name = d + "-" + a;
          final Integer score = mFusionScores.get(name);
          if (score != null && score > bestScore) {
            bestFusion = name;
            bestScore = score;
          }
        }
      }
      if (bestScore > 0) {
        rec.setInfo(INFO_FUSION_NAME, bestFusion);
        rec.setInfo(INFO_FUSION_SCORE, bestScore.toString());
      }
    }
  }

  // Drops records that could not produce a transcript between an acceptor and donor gene
  private static class RegionFilter implements VcfFilter {
    Pair<ReferenceRanges<String>, ReferenceRanges<String>> mAcceptorRegions;
    Pair<ReferenceRanges<String>, ReferenceRanges<String>> mDonorRegions;
    RegionFilter(File acceptorBed, File donorBed, int geneCol, int strandCol) throws IOException {
      this(loadGeneBed(acceptorBed, geneCol, strandCol), loadGeneBed(donorBed, geneCol, strandCol));
    }
    RegionFilter(Pair<ReferenceRanges<String>, ReferenceRanges<String>> acceptorRegions,
                 Pair<ReferenceRanges<String>, ReferenceRanges<String>> donorRegions) {
      mAcceptorRegions = acceptorRegions;
      mDonorRegions = donorRegions;
    }

    @Override
    public boolean accept(VcfRecord record) {
      final VariantType vt = VariantType.getType(record.getAllele(0), record.getAllele(1));
      switch (vt) {
        case SV_BREAKEND:
          return possibleTranscript(record.getSequenceName(), record.getStart(), new BreakpointAlt(record.getAltCalls().get(0)), record);
        case SV_SYMBOLIC:
          final ArrayList<String> svTypes = record.getInfo().get(INFO_SVTYPE);
          if (svTypes == null || svTypes.size() != 1) {
            return false;
          }
          switch (svTypes.get(0)) {
            case "DEL":
              // Make equivalent breakend for the non-deleted section. Only one breakend required
              return possibleTranscript(record.getSequenceName(), record.getStart(), new BreakpointAlt(record.getRefCall() + "[" + record.getSequenceName() + ":" + VcfUtils.getIntegerInfoFieldFromRecord(record, "END") + "["), record);

            case "INV":
              // Make equivalent breakends for both continuations of the current position. Far end not needed as it's checked internally
              return possibleTranscript(record.getSequenceName(), record.getStart(), new BreakpointAlt(record.getRefCall() + "]" + record.getSequenceName() + ":" + VcfUtils.getIntegerInfoFieldFromRecord(record, "END") + "]"), record)
                || possibleTranscript(record.getSequenceName(), record.getStart() + 1, new BreakpointAlt("[" + record.getSequenceName() + ":" + (VcfUtils.getIntegerInfoFieldFromRecord(record, "END") + 1) + "[" + record.getRefCall()), record);
            default:
            // DUP cannot result in a fusion
            // INS are typically small-scale and also can't yield a fusion
              return false;
          }
        default:
          break;
      }
      return false;
    }

    /**
     * Determines whether the breakpoint could represent a transcript between donor and acceptor genes.
     * Checks both strands.
     * @param chr chromosome of the break point
     * @param pos position of the break point
     * @param breakend the breakend encoding remote position and orientations
     * @param record to be annotated with the name of acceptor and donor genes, if any
     * @return true if the breakpoint represents a matching transcript.
     */
    private boolean possibleTranscript(String chr, int pos, BreakpointAlt breakend, VcfRecord record) {
      // Check for a transcript in either direction
      return possibleTranscript(record,
        getGenes(breakend.isLocalUp() ? mAcceptorRegions.getB() : mAcceptorRegions.getA(), chr, pos),
        getGenes(breakend.isRemoteUp() ? mDonorRegions.getA() : mDonorRegions.getB(), breakend.getRemoteChr(), breakend.getRemotePos()))
        || possibleTranscript(record,
        getGenes(breakend.isRemoteUp() ? mAcceptorRegions.getB() : mAcceptorRegions.getA(), breakend.getRemoteChr(), breakend.getRemotePos()),
        getGenes(breakend.isLocalUp() ? mDonorRegions.getA() : mDonorRegions.getB(), chr, pos));
    }

    private boolean possibleTranscript(VcfRecord record, List<String> acceptorGenes, List<String> donorGenes) {
      if (acceptorGenes != null && donorGenes != null) {
        record.setInfo(INFO_ACCEPTOR_GENE, acceptorGenes.toArray(new String[acceptorGenes.size()]));
        record.setInfo(INFO_DONOR_GENE, donorGenes.toArray(new String[donorGenes.size()]));
        return true;
      }
      return false;
    }

    private List<String> getGenes(ReferenceRanges<String> ranges, String sequence, int position) {
      final RangeList<String> l = ranges.get(sequence);
      if (l != null) {
        return l.find(position);
      }
      return null;
    }

    @Override
    public void setHeader(VcfHeader header) {
      header.ensureContains(new InfoField(INFO_ACCEPTOR_GENE, MetaType.STRING, VcfNumber.DOT, "Name of acceptor gene(s)"));
      header.ensureContains(new InfoField(INFO_DONOR_GENE, MetaType.STRING, VcfNumber.DOT, "Name of donor gene(s)"));
    }

    private static Pair<ReferenceRanges<String>, ReferenceRanges<String>> loadGeneBed(File bedFile, int geneCol, int strandCol) throws IOException {
      final ReferenceRanges.Accumulator<String> posRangesAccum = new ReferenceRanges.Accumulator<>();
      final ReferenceRanges.Accumulator<String> negRangesAccum = new ReferenceRanges.Accumulator<>();
      try (BedReader bedReader = BedReader.openBedReader(null, bedFile, Math.max(geneCol, strandCol) + 1)) {
        while (bedReader.hasNext()) {
          final BedRecord record = bedReader.next();
          final String gene = record.getAnnotations()[geneCol];
          final String strand = record.getAnnotations()[strandCol];
          final RangeList.RangeData<String> rangeData = new RangeList.RangeData<>(record.getStart(), record.getEnd(), gene);
          if ("+".equals(strand)) {
            posRangesAccum.addRangeData(record.getSequenceName(), rangeData);
          } else if ("-".equals(strand)) {
            negRangesAccum.addRangeData(record.getSequenceName(), rangeData);
          } else {
            throw new IOException("BED record does not indicate strand as '-' or '+' in column 6: " + record);
          }
        }
      }
      return new Pair<>(posRangesAccum.getReferenceRanges(), negRangesAccum.getReferenceRanges());
    }
  }
}
