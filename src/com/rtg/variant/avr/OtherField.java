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
package com.rtg.variant.avr;

import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.TypedField;
import com.rtg.vcf.header.VcfNumber;

/**
 * Dummy TypedField for AVR fields that do not come from INFO or FORMAT fields
 */
public class OtherField implements TypedField<OtherField> {

  /** Singleton instance there need only be one */
  public static final OtherField SINGLETON = new OtherField();

  @Override
  public VcfNumber getNumber() {
    return VcfNumber.ONE;
  }

  @Override
  public MetaType getType() {
    return MetaType.FLOAT;
  }

  @Override
  public String getDescription() {
    return "OTHER";
  }

  @Override
  public String getId() {
    return "OTHER";
  }

  @Override
  public OtherField superSet(OtherField other) {
    return this;
  }
}
