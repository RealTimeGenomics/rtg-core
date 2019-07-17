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

package com.rtg.variant.sv.discord;

import static com.rtg.vcf.VcfUtils.FILTER_PASS;
import static com.rtg.vcf.VcfUtils.FORMAT_GENOTYPE;
import static com.rtg.vcf.VcfUtils.INFO_CIPOS;
import static com.rtg.vcf.VcfUtils.INFO_COMBINED_DEPTH;
import static com.rtg.vcf.VcfUtils.INFO_IMPRECISE;
import static com.rtg.vcf.VcfUtils.INFO_SVTYPE;
import static com.rtg.vcf.VcfUtils.SvType.BND;

import java.io.IOException;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

import com.rtg.mode.DNA;
import com.rtg.reader.NamesInterface;
import com.rtg.reader.SequencesReader;
import com.rtg.vcf.BreakpointAlt;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

import htsjdk.samtools.SAMFileHeader;

/**
 * Assumes a calling pattern that presents breakends for each template in groups.
 */
public class VcfDiscordantOutputFormatter {

  private static final String INFO_COVERAGE = "CV";
  private static final String INFO_AMBIGUITY = "AR";
  private static final String FILTER_INCONSISTENT = "INCONSISTENT";

  private final SequencesReader mTemplate;

  private final Map<String, Integer> mSequenceMap = new HashMap<>();

  private String mCurrentSequenceName;
  private byte[] mCurrentSequence;


  /**
   * Create a new object
   * @param genomeSequences reader for template
   * @throws IOException if error
   */
  public VcfDiscordantOutputFormatter(SequencesReader genomeSequences) throws IOException {
    mTemplate = genomeSequences;
    final NamesInterface pni = genomeSequences.names();
    for (long i = 0; i < pni.length(); ++i) {
      mSequenceMap.put(genomeSequences.names().name(i), (int) i);
    }
  }

  /**
   * Given a breakpoint constraint, return VCF
   *
   * @param readset discordant read set
   * @param coverage the coverage at this position, -1 if you'd like to avoid output
   * @param ambiguous the ration between ambiguous and unambiguous reads, -1 if you'd like to avoid output
   * @return {@link VcfRecord} containing current results
   *
   * @throws IOException if error
   */
  public VcfRecord vcfRecord(DiscordantReadSet readset, int coverage, double ambiguous) throws IOException {
    final boolean unionOnly = readset.getIntersection() == null;
    final BreakpointConstraint geo = unionOnly ? readset.getUnion() : readset.getIntersection();
    final BreakpointPosition pos = geo.position();
    final int ory = geo.getOrientation().yDir();
    final int orx = geo.getOrientation().xDir();
    // Adjust ref pos include anchoring base for "local up" breakends
    // Use Math.max(x, 0) to handle breakends positioned before the start of the reference
    final int refAdjust = orx == +1 ? 1 : 0;
    final int refPosition = Math.max(pos.position() - refAdjust, 1 - refAdjust);
    final String cipos = "" + (pos.lo() - refPosition - refAdjust) + "," + (pos.hi() - refPosition - refAdjust);
    final String ref = getRef(geo.getXName(), refPosition);
    final VcfRecord rec = new VcfRecord(geo.getXName(), refPosition, ref);
    rec.setId(VcfRecord.MISSING);
    final String alt = new BreakpointAlt(ref, orx == +1, geo.getYName(), Math.max(pos.positionAlt(), 0), ory == +1).toString();
    rec.addAltCall(alt);
    rec.setQuality(VcfRecord.MISSING);
    if (unionOnly) {
      rec.addFilter(FILTER_INCONSISTENT);
    } else {
      rec.addFilter(FILTER_PASS);
    }
    rec.setInfo(INFO_IMPRECISE);
    rec.setInfo(INFO_SVTYPE, BND.name());
    rec.setInfo(INFO_COMBINED_DEPTH, "" + readset.getCounts());
    rec.setInfo(INFO_CIPOS, cipos);
    if (coverage > -1) {
      rec.setInfo(INFO_COVERAGE, "" + coverage);
    }
    if (ambiguous > -1) {
      rec.setInfo(INFO_AMBIGUITY, "" + ambiguous);
    }
    rec.setNumberOfSamples(1);
    rec.addFormatAndSample(FORMAT_GENOTYPE, "1/1"); //for now default to homozygous
    return rec;
  }

  /**
   * @param refName reference name
   * @param pos one based position
   * @return nt at asked position
   */
  private String getRef(String refName, int pos) throws IOException {
    if (!refName.equals(mCurrentSequenceName)) {
      final Integer id = mSequenceMap.get(refName);
      if (id == null) {
        throw new RuntimeException("reference name was not contained in given SDF : " + refName);
      }
      mCurrentSequenceName = refName;
      mCurrentSequence = new byte[mTemplate.length(id)];
      mTemplate.read(id, mCurrentSequence);
    }
    assert pos < mCurrentSequence.length : "pos=" + pos + " currentLen=" + mCurrentSequence.length;
    final DNA dna = (pos < 0 || (pos > mCurrentSequence.length - 1)) ? DNA.N : DNA.values()[mCurrentSequence[pos]];
    return dna.toString().toUpperCase(Locale.getDefault());
  }

  /**
   * Returns the VcfHeader
   * @param samheader a SAM header, may be null
   * @param samplename the sample name
   * @param coverage include coverage header
   * @param ambiguity include ambiguity header
   * @return the VcfHeader
   */
  public VcfHeader header(SAMFileHeader samheader, String samplename, boolean coverage, boolean ambiguity) {
    final VcfHeader header = new VcfHeader();
    header.addCommonHeader();
    header.addReference(mTemplate);
    if (samheader != null) {
      header.addContigFields(samheader);
    }
    header.addFilterField(FILTER_INCONSISTENT, "Supporting reads are inconsistent as to breakend location");
    header.addInfoField(INFO_CIPOS, MetaType.INTEGER, new VcfNumber("2"), "Confidence interval around POS for imprecise variants");
    header.addInfoField(INFO_IMPRECISE, MetaType.FLAG, new VcfNumber("0"), "Imprecise structural variation");
    header.addInfoField(INFO_SVTYPE, MetaType.STRING, VcfNumber.ONE, "Type of structural variant");
    header.addInfoField(INFO_COMBINED_DEPTH, MetaType.INTEGER, VcfNumber.ONE, "Read Depth");
    if (coverage) {
      header.addInfoField(INFO_COVERAGE, MetaType.INTEGER, VcfNumber.ONE, "Coverage at start position");
    }
    if (ambiguity) {
      header.addInfoField(INFO_AMBIGUITY, MetaType.INTEGER, VcfNumber.ONE, "Ambiguity at start position");
    }
    header.addFormatField(FORMAT_GENOTYPE, MetaType.STRING, VcfNumber.ONE, "Genotype");
    header.addSampleName(samplename);
    return header;
  }

}
