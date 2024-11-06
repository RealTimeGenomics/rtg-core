/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package com.rtg.vcf.validator;

import com.rtg.vcf.VcfReaderTest;
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
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=0.0\tGT\t0/1"));
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=1.243\tGT\t0/1"));
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=2e10\tGT\t0/1"));
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=3\tGT\t0/1"));
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=.\tGT\t0/1"));
    } catch (RuleValidationException e) {
      fail("These records should be passing: " + e.getMessage());
    }
    try {
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=bar\tGT\t0/1"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("Value in INFO field FOO not in the correct number format.", e.getMessage());
    }
    try {
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=NaN\tGT\t0/1"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("One or more values for the INFO field FOO is outside the expected range of values.", e.getMessage());
    }
    try {
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=Infinity\tGT\t0/1"));
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
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=2\tGT\t0/1"));
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=10\tGT\t0/1"));
    } catch (RuleValidationException e) {
      fail("These records should be passing: " + e.getMessage());
    }
    try {
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=1\tGT\t0/1"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("One or more values for the INFO field FOO is outside the expected range of values.", e.getMessage());
    }
    try {
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=0\tGT\t0/1"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("One or more values for the INFO field FOO is outside the expected range of values.", e.getMessage());
    }
    try {
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=20\tGT\t0/1"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("One or more values for the INFO field FOO is outside the expected range of values.", e.getMessage());
    }
    try {
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=21\tGT\t0/1"));
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
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=2\tGT\t0/1"));
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=10\tGT\t0/1"));
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=1\tGT\t0/1"));
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=20\tGT\t0/1"));
    } catch (RuleValidationException e) {
      fail("These records should be passing: " + e.getMessage());
    }
    try {
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=0\tGT\t0/1"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("One or more values for the INFO field FOO is outside the expected range of values.", e.getMessage());
    }
    try {
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=21\tGT\t0/1"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("One or more values for the INFO field FOO is outside the expected range of values.", e.getMessage());
    }
  }
}
