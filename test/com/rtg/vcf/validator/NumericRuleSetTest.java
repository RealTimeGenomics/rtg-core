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

package com.rtg.vcf.validator;

import com.rtg.vcf.VcfReader;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfNumber;
import com.rtg.vcf.validator.NumericRuleSet.DoubleConverter;
import com.rtg.vcf.validator.NumericRuleSet.LongConverter;
import com.rtg.vcf.validator.RuleSet.FieldType;

import junit.framework.TestCase;

/**
 */
public class NumericRuleSetTest extends TestCase {

  public void testLongConverter() {
    final LongConverter conv = new LongConverter();
    assertEquals(Long.valueOf(8), conv.getValue("8"));
  }

  public void testDoubleConverter() {
    final DoubleConverter conv = new DoubleConverter();
    assertEquals(8.56, conv.getValue("8.56"));
  }

  public void testNumericConversionValidation() {
    final NumericRuleSet<Double> set = new NumericRuleSet<>("FOO", FieldType.INFO, VcfNumber.ONE, MetaType.FLOAT, new DoubleConverter());
    set.addNaNRule();
    set.addInfinityRule();
    try {
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=0.0\tGT\t0/1"));
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=1.243\tGT\t0/1"));
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=2e10\tGT\t0/1"));
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=3\tGT\t0/1"));
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=.\tGT\t0/1"));
    } catch (RuleValidationException e) {
      fail("These records should be passing: " + e.getMessage());
    }
    try {
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=bar\tGT\t0/1"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("Value in INFO field FOO not in the correct number format.", e.getMessage());
    }
    try {
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=NaN\tGT\t0/1"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("One or more values for the INFO field FOO is outside the expected range of values.", e.getMessage());
    }
    try {
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=Infinity\tGT\t0/1"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("One or more values for the INFO field FOO is outside the expected range of values.", e.getMessage());
    }
  }

  public void testLessGreaterThanRules() throws RuleValidationException {
    final NumericRuleSet<Long> set = new NumericRuleSet<>("FOO", FieldType.INFO, VcfNumber.ONE, MetaType.INTEGER, new LongConverter());
    set.addGreaterThanRule("1");
    set.addLessThanRule("20");
    try {
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=2\tGT\t0/1"));
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=10\tGT\t0/1"));
    } catch (RuleValidationException e) {
      fail("These records should be passing: " + e.getMessage());
    }
    try {
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=1\tGT\t0/1"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("One or more values for the INFO field FOO is outside the expected range of values.", e.getMessage());
    }
    try {
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=0\tGT\t0/1"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("One or more values for the INFO field FOO is outside the expected range of values.", e.getMessage());
    }
    try {
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=20\tGT\t0/1"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("One or more values for the INFO field FOO is outside the expected range of values.", e.getMessage());
    }
    try {
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=21\tGT\t0/1"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("One or more values for the INFO field FOO is outside the expected range of values.", e.getMessage());
    }
  }

  public void testLessGreaterThanEqualRules() throws RuleValidationException {
    final NumericRuleSet<Long> set = new NumericRuleSet<>("FOO", FieldType.INFO, VcfNumber.ONE, MetaType.INTEGER, new LongConverter());
    set.addGreaterThanEqualRule("1");
    set.addLessThanEqualRule("20");
    try {
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=2\tGT\t0/1"));
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=10\tGT\t0/1"));
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=1\tGT\t0/1"));
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=20\tGT\t0/1"));
    } catch (RuleValidationException e) {
      fail("These records should be passing: " + e.getMessage());
    }
    try {
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=0\tGT\t0/1"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("One or more values for the INFO field FOO is outside the expected range of values.", e.getMessage());
    }
    try {
      set.validateRecord(VcfReader.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=21\tGT\t0/1"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("One or more values for the INFO field FOO is outside the expected range of values.", e.getMessage());
    }
  }
}
