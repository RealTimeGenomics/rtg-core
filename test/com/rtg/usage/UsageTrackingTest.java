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
package com.rtg.usage;

import java.io.File;
import java.io.IOException;
import java.util.Date;
import java.util.Properties;
import java.util.UUID;

import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;

import junit.framework.TestCase;

/**
 */
public class UsageTrackingTest extends TestCase {
  public void testUsageTracking() throws IOException {
    final Properties prop = new Properties();
    prop.setProperty(UsageTracking.REQUIRE_USAGE, "false");
    prop.setProperty(UsageTracking.USAGE_DESTINATION, UsageTracking.UsageDestination.FILE_OR_SERVER.toString());
    try (final TestDirectory dir = new TestDirectory()) {
      final File cfgFile = UsageConfiguration.createSimpleConfigurationFile(new File(dir, "cfgFile"), dir.getPath(), null);
      final UUID runId = UUID.randomUUID();
      final UsageTracking usage = new UsageTracking(prop, "testModule", runId, cfgFile, false);
      usage.recordBeginning();
      usage.recordEnd(428428791, true);
      final File usageFile = FileUsageTrackingClient.ensureUsageFile(dir, new Date());
      final String res = FileUtils.fileToString(usageFile);
      TestUtils.containsAll(res, "testModule", "428428791", runId.toString(), "Success", "Start");
      assertTrue(usage.allowFileTracking());
      assertFalse(usage.requireUsage());
    }
  }
}
