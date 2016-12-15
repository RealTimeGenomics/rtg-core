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

import java.util.Collections;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.vcf.VcfRecord;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class CnvRecordFilterTest extends TestCase {

  public void test() {
    Diagnostic.setLogStream();
    final CnvRecordFilter f = new CnvRecordFilter(Collections.singletonList("pretend"));
    final VcfRecord record = new VcfRecord("pretend", 42, "A");
    assertFalse(f.accept(record));
    record.addInfo(CnaType.INFO_END, "42");
    assertFalse(f.accept(record));
    record.addInfo(CnaType.INFO_SVTYPE, CnaType.DEL.toString());
    assertTrue(f.accept(record));
    final VcfRecord record2 = new VcfRecord("pretend2", 42, "A");
    record2.addInfo(CnaType.INFO_END, "42");
    record2.addInfo(CnaType.INFO_SVTYPE, CnaType.DEL.toString());
    assertFalse(f.accept(record2));
  }

}
