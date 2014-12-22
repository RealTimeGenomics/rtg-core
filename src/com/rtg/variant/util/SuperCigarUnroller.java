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
package com.rtg.variant.util;

import com.rtg.mode.DNA;
import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.SuperCigarParser;
import com.rtg.variant.VariantAlignmentRecord;

/**
 * Validates a super cigar against the template and original read.
 * This reconstructs the original read from the alignment record plus the template.
 * Usage example:
 * <pre>
 *   SuperCigarUnroller unroll = new SuperCigarUnroller();
 *   unroll.setSamRecord(samrec);
 *   unroll.setTemplate(template);
 *   unroll.parse();
 *   String originalRead = unroll.getString();
 * </pre>
 *
 */
public class SuperCigarUnroller extends SuperCigarParser {

  protected StringBuilder mBuilder;

  /**
   * Set the SAM record, which must contain a super cigar and read delta attribute.
   * @param alignmentRecord a valid SAM record.
   */
  public void setAlignmentRecord(VariantAlignmentRecord alignmentRecord) {
    final String superCigar = alignmentRecord.getSuperCigar();
    final String cigar =  superCigar == null ? alignmentRecord.getCigar() : superCigar;
    super.setCigar(cigar, alignmentRecord.getCGReadDelta());
    setTemplateStart(alignmentRecord.getStart());
    mBuilder = new StringBuilder();
  }


  /**
   * @return the read string, with spaces for CG gaps.
   */
  public String getString() {
    return mBuilder.toString();
  }

  /**
   * @throws BadSuperCigarException on bad cigar.
   */
  @Override
  public void parse() throws BadSuperCigarException {
    mBuilder = new StringBuilder();
    super.parse();
    //checkMismatchCounts();
  }

  @Override
  protected void doReadHardClip() {
  }

  @Override
  protected void doReadSoftClip(int readNt) {
    mBuilder.append(DNA.valueChars()[readNt]);
  }

  @Override
  protected void doTemplateOverlap() {
  }

  @Override
  protected void doTemplateSkip(int templateNt) {
  }

  @Override
  protected void doReadOnly(int readNt) {
    mBuilder.append(DNA.valueChars()[readNt]);
  }

  @Override
  protected void doTemplateOnly(int templateNt) {
  }

  @Override
  protected void doSubstitution(int readNt, int templateNt) {
    mBuilder.append(DNA.valueChars()[readNt]);
  }

  @Override
  protected void doEquality(int readNt, int nt) {
    mBuilder.append(DNA.valueChars()[nt]);
  }


  @Override
  protected void doUnknownOnTemplate(int readNt, int templateNt) {
    mBuilder.append(DNA.valueChars()[readNt]);
  }


  @Override
  protected void doUnknownOnRead() {
    mBuilder.append(DNA.N.name());
  }


}
