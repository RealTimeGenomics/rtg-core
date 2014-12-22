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

import java.io.IOException;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

import com.rtg.mode.DNA;
import com.rtg.reader.PrereadNamesInterface;
import com.rtg.reader.SequencesReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

import net.sf.samtools.SAMFileHeader;

/**
 * Assumes a calling pattern that presents break ends for each template in groups.
 */
public class VcfDiscordantOutputFormatter {

  //private static final String INFO_DIRECTION = "DIRECTION";
  private static final String INFO_CIPOS = "CIPOS";
  private static final String INFO_DP = "DP";
  private static final String INFO_COVERAGE = "CV";
  private static final String INFO_AMBIGUITY = "AR";
  private static final String INFO_SVTYPE = "SVTYPE";
  private static final String INFO_IMPRECISE = "IMPRECISE";
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
    final PrereadNamesInterface pni = genomeSequences.names();
    for (long i = 0; i < pni.length(); i++) {
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
    final VcfRecord rec = new VcfRecord();
    rec.setSequence(geo.getXName());
    final BreakpointPosition pos = geo.position();
    //Use Math.max(x, 0) to fix records that go before the start of the reference
    final int refPosition = Math.max(pos.position(), 0);
    final String cipos = "" + (pos.lo() - refPosition) + "," + (pos.hi() - refPosition);
    rec.setStart(refPosition);
    rec.setId(VcfRecord.MISSING);
    final String ref = getRef(geo.getXName(), refPosition);
    rec.setRefCall(ref);
    final String alt;

    final String bracket;
    final int ory = geo.getOrientation().getY();
    if (ory == +1) {
      bracket = "]";
    } else {
      assert ory == -1;
      bracket = "[";
    }
    final String alt0 = bracket + geo.getYName() + ":" + (Math.max(pos.positionAlt(), 0) + 1) + bracket;
    final int orx = geo.getOrientation().getX();
    if (orx == +1) {
      alt = ref + alt0;
    } else {
      assert orx == -1;
      alt = alt0 + ref;
    }

    rec.addAltCall(alt);
    rec.setQuality(VcfRecord.MISSING);
    if (unionOnly) {
      rec.addFilter(FILTER_INCONSISTENT);
    } else {
      rec.addFilter(VcfUtils.FILTER_PASS);
    }
    rec.addInfo(INFO_IMPRECISE);
    rec.addInfo(INFO_SVTYPE, "BND");
    rec.addInfo(INFO_DP, "" + readset.getCounts());
    rec.addInfo(INFO_CIPOS, cipos);
    if (coverage > -1) {
      rec.addInfo(INFO_COVERAGE, "" + coverage);
    }
    if (ambiguous > -1) {
      rec.addInfo(INFO_AMBIGUITY, "" + ambiguous);
    }
    //rec.addInfo(INFO_DIRECTION, geo.getOrientation().name());
    rec.setNumberOfSamples(1);
    rec.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, "1/1"); //for now default to homozygous
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
      mTemplate.seek(id);
      mCurrentSequence = new byte[mTemplate.currentLength()];
      mTemplate.readCurrent(mCurrentSequence);
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
    //header.addInfoField(INFO_DIRECTION, MetaType.STRING, new VcfNumber("1"), "Direction for current breakpoint");
    header.addInfoField(INFO_SVTYPE, MetaType.STRING, new VcfNumber("1"), "Type of structural variant");
    header.addInfoField(INFO_DP, MetaType.INTEGER, new VcfNumber("1"), "Read Depth");
    if (coverage) {
      header.addInfoField(INFO_COVERAGE, MetaType.INTEGER, new VcfNumber("1"), "Coverage at start position");
    }
    if (ambiguity) {
      header.addInfoField(INFO_AMBIGUITY, MetaType.INTEGER, new VcfNumber("1"), "Ambiguity at start position");
    }
    header.addFormatField(VcfUtils.FORMAT_GENOTYPE, MetaType.STRING, new VcfNumber("1"), "Genotype");
    header.addSampleName(samplename);
    return header;
  }

}
