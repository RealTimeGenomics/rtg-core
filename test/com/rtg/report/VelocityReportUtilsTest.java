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
package com.rtg.report;

import java.io.IOException;
import java.util.HashMap;

import com.rtg.util.HtmlReportHelper;
import com.rtg.util.TestUtils;
import com.rtg.util.io.TestDirectory;

import junit.framework.TestCase;

/**
 */
public class VelocityReportUtilsTest extends TestCase {

  public void testTemplate() throws IOException {
    final HashMap<String, String> map = new HashMap<>();
    map.put("aVariable", "foorah");
    map.put("bVariable", "inner stuff");
    final String res = VelocityReportUtils.processTemplate("velocityTest.txt", map);
    assertEquals("foorah\ntextNextToinner stufffoo\n\n", res);
    map.put("valExists", "all new value");
    final String res2 = VelocityReportUtils.processTemplate("velocityTest.txt", map);
    assertEquals("foorah\ntextNextToinner stufffoo\n\nall new value\n", res2);
  }

  public void testWrapDefault() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final HtmlReportHelper hrh = new HtmlReportHelper(dir, "index");
      final String res = VelocityReportUtils.wrapDefaultTemplate("hardbody", "test", hrh);

      TestUtils.containsAll(res, "\nhardbody\n    </div>",
                                   "title>test</title",
                                   "href=\"" + hrh.getResourcesDirName() + "/rtg.css");
    }
  }
}
