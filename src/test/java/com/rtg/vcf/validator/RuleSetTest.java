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

import com.rtg.util.TestUtils;
import com.rtg.vcf.VcfReaderTest;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfNumber;
import com.rtg.vcf.header.VcfNumberType;
import com.rtg.vcf.validator.RuleSet.FieldType;
import com.rtg.vcf.validator.RuleSet.StringConverter;

import junit.framework.TestCase;

/**
 */
public class RuleSetTest extends TestCase {

  public void testFieldType() {
    TestUtils.testEnum(FieldType.class, "[INFO, FORMAT]");
  }

  public void testStringConverter() {
    final String str = "testString";
    final StringConverter conv = new StringConverter();
    assertEquals(str, conv.getValue(str));
  }

  public void testEnumerationRule() {
    final RuleSet<String> set = new RuleSet<>("FOO", FieldType.INFO, VcfNumber.ONE, MetaType.STRING, new StringConverter());
    try {
      set.addEnumerationRule("a,b,c");
    } catch (RuleValidationException e) {
      fail("Should not be able to fail a string to string conversion");
    }
    try {
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=a\tGT\t0/1"));
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=b\tGT\t0/1"));
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=c\tGT\t0/1"));
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=.\tGT\t0/1"));
    } catch (RuleValidationException e) {
      fail("These records should be passing: " + e.getMessage());
    }
    try {
      set.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=bar\tGT\t0/1"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("One or more values for the INFO field FOO is outside the expected range of values.", e.getMessage());
    }
  }

  public void testGetters() {
    final RuleSet<String> set = new RuleSet<>("FOO", FieldType.FORMAT, new VcfNumber("2"), MetaType.STRING, new StringConverter());
    assertEquals("FOO", set.getName());
    assertEquals(FieldType.FORMAT, set.getVcfFieldType());
    assertEquals(2, set.getVcfNumber().getNumber());
    assertEquals(VcfNumberType.INTEGER, set.getVcfNumber().getNumberType());
    assertEquals(MetaType.STRING, set.getMetaType());
  }

  public void testMultipleValueValidation() {
    final RuleSet<String> fooSet = new RuleSet<>("FOO", FieldType.INFO, new VcfNumber("A"), MetaType.STRING, new StringConverter());
    try {
      fooSet.addEnumerationRule("a,b,c");
    } catch (RuleValidationException e) {
      fail("Should not be able to fail a string to string conversion");
    }
    try {
      fooSet.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=a\tGT\t0/1"));
      fooSet.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT,G\t167.9\tPASS\tFOO=a,c\tGT\t1/2"));
      fooSet.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT,G\t167.9\tPASS\tFOO=.,c\tGT\t1/2"));
      fooSet.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\t.\t167.9\tPASS\tFOO=.\tGT\t1/2"));
    } catch (RuleValidationException e) {
      fail("These records should be passing: " + e.getMessage());
    }
    try {
      fooSet.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT,G\t167.9\tPASS\tFOO=a,bar\tGT\t1/2"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("One or more values for the INFO field FOO is outside the expected range of values.", e.getMessage());
    }
    try {
      fooSet.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT,G\t167.9\tPASS\tFOO=a,b,c\tGT\t1/2"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("Too many values for the INFO field FOO.", e.getMessage());
    }
    try {
      fooSet.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT,G\t167.9\tPASS\tFOO=a\tGT\t1/2"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("Too few values for the INFO field FOO.", e.getMessage());
    }
    try {
      fooSet.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT,G\t167.9\tPASS\tFOO\tGT\t1/2"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("Expected values for the INFO field FOO.", e.getMessage());
    }
  }

  public void testAnyNumberValidation() {
    final RuleSet<String> fooSet = new RuleSet<>("FOO", FieldType.INFO, new VcfNumber("."), MetaType.STRING, new StringConverter());
    try {
      fooSet.addEnumerationRule("a,b,c");
    } catch (RuleValidationException e) {
      fail("Should not be able to fail a string to string conversion");
    }
    try {
      fooSet.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=a\tGT\t0/1"));
      fooSet.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=a,c,b\tGT\t1/2"));
      fooSet.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=.,c\tGT\t1/2"));
      fooSet.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\t.\t167.9\tPASS\tFOO=.\tGT\t1/2"));
    } catch (RuleValidationException e) {
      fail("These records should be passing: " + e.getMessage());
    }
  }

  public void testFlagType() {
    final RuleSet<String> fooSet = new RuleSet<>("FOO", FieldType.INFO, VcfNumber.FLAG, MetaType.FLAG, null);
    try {
      fooSet.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO\tGT\t0/1"));
    } catch (RuleValidationException e) {
      fail("These records should be passing: " + e.getMessage());
    }
    try {
      fooSet.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\tFOO=a\tGT\t0/1"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("Expected no values for the INFO field FOO.", e.getMessage());
    }
  }

  public void testFormatValidation() {
    final RuleSet<String> fooSet = new RuleSet<>("FOO", FieldType.FORMAT, new VcfNumber("2"), MetaType.STRING, new StringConverter());
    try {
      fooSet.addEnumerationRule("a,b,c");
    } catch (RuleValidationException e) {
      fail("Should not be able to fail a string to string conversion");
    }
    try {
      fooSet.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\t.\tGT:FOO\t0/1:a,b\t0/0:c,.\t1/1:."));
    } catch (RuleValidationException e) {
      fail("These records should be passing: " + e.getMessage());
    }
    try {
      fooSet.validateRecord(VcfReaderTest.vcfLineToRecord("chr21\t1024\t.\tC\tT\t167.9\tPASS\t.\tGT:FOO\t0/1:a,b\t0/0:c,.\t1/1:bar,b"));
      fail("This record should have failed.");
    } catch (RuleValidationException e) {
      assertEquals("One or more values for the FORMAT field FOO is outside the expected range of values.", e.getMessage());
    }
  }

}
