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
import com.rtg.vcf.VcfUtils;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class CnvRecordFilterTest extends TestCase {

  public void test() {
    Diagnostic.setLogStream();
    final CnvRecordFilter f = new CnvRecordFilter(Collections.singletonList("pretend"), true);
    f.setHeader(null);
    final VcfRecord record = new VcfRecord("pretend", 42, "A");
    assertFalse(f.accept(record)); // Not SV
    record.addInfo(VcfUtils.INFO_END, "42");
    assertFalse(f.accept(record)); // Not SV
    record.addInfo(VcfUtils.INFO_SVTYPE, CnaType.DEL.toString());
    assertTrue(f.accept(record)); // SV with END

    final VcfRecord record2 = new VcfRecord("pretend2", 42, "A");
    record2.addInfo(VcfUtils.INFO_END, "42");
    record2.addInfo(VcfUtils.INFO_SVTYPE, CnaType.DEL.toString());
    assertFalse(f.accept(record2)); // Overlap
  }

  public void testNoEnd() {
    Diagnostic.setLogStream();
    final CnvRecordFilter f = new CnvRecordFilter(Collections.singletonList("pretend"), true);
    final VcfRecord record = new VcfRecord("pretend", 42, "A");
    record.addInfo(VcfUtils.INFO_SVTYPE, CnaType.DEL.toString());
    assertFalse(f.accept(record));
  }

  public void testOverlap() {
    Diagnostic.setLogStream();
    final CnvRecordFilter f = new CnvRecordFilter(Collections.singletonList("pretend"), true);
    final VcfRecord record = new VcfRecord("pretend", 42, "A");
    record.addInfo(VcfUtils.INFO_END, "48");
    record.addInfo(VcfUtils.INFO_SVTYPE, CnaType.DEL.toString());
    assertTrue(f.accept(record));

    final VcfRecord record2 = new VcfRecord("pretend", 46, "A");
    record2.addInfo(VcfUtils.INFO_END, "49");
    record2.addInfo(VcfUtils.INFO_SVTYPE, CnaType.DEL.toString());
    assertFalse(f.accept(record2));

    final VcfRecord record3 = new VcfRecord("pretend", 47, "A"); // Since according to VCF spec, SV records include the base BEFORE the event
    record3.addInfo(VcfUtils.INFO_END, "49");
    record3.addInfo(VcfUtils.INFO_SVTYPE, CnaType.DEL.toString());
    assertTrue(f.accept(record3));
  }

}
