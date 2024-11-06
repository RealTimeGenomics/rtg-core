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

package com.rtg.variant.avr;

import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.annotation.AbstractDerivedInfoAnnotation;
import com.rtg.vcf.header.InfoField;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

import junit.framework.TestCase;

/**
 */
public class DerivedAnnotationTest extends TestCase {

  public void testConstructor() {
    final DerivedAnnotation ann = new DerivedAnnotation("IC");
    assertEquals("DERIVED-IC", ann.getName());
    assertEquals(AnnotationDataType.DOUBLE, ann.getType());
  }

  private static final class DummyVcfAnnotation extends AbstractDerivedInfoAnnotation {
    DummyVcfAnnotation() {
      super(new InfoField("DUMMY", MetaType.FLAG, VcfNumber.ONE, "DUMMY-DESC"), null);
    }
    @Override
    public Object getValue(VcfRecord record, int sampleNumber) {
      assertNull(record);
      return true;
    }
    @Override
    public String checkHeader(VcfHeader header) {
      assertNull(header);
      return "A String";
    }
  }

  public void testWrapping() {
    final DerivedAnnotation ann = new DerivedAnnotation(new DummyVcfAnnotation());
    assertEquals("DERIVED-DUMMY", ann.getName());
    assertEquals(AnnotationDataType.BOOLEAN, ann.getType());
    assertEquals(true, ann.getValue(null, 0));
    assertEquals("A String", ann.checkHeader(null));
  }
}
