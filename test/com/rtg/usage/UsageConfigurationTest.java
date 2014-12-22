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
import java.util.Properties;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.TestDirectory;

import junit.framework.TestCase;

/**
 */
public class UsageConfigurationTest extends TestCase {
  static UsageConfiguration createSimpleConfiguration(String usageDir, String usageHost, Boolean logUsername, Boolean logHostname, Boolean logCommandline) {
    final Properties prop = new Properties();
    if (usageDir != null) {
      prop.setProperty(UsageConfiguration.USAGE_DIR, usageDir);
    }
    if (usageHost != null) {
      prop.setProperty(UsageConfiguration.USAGE_HOST, usageHost);
    }
    if (logUsername != null) {
      prop.setProperty(UsageConfiguration.USAGE_LOG_USERNAME, logUsername.toString());
    }
    if (logHostname != null) {
      prop.setProperty(UsageConfiguration.USAGE_LOG_HOSTNAME, logHostname.toString());
    }
    if (logCommandline != null) {
      prop.setProperty(UsageConfiguration.USAGE_LOG_COMMANDLINE, logCommandline.toString());
    }
    return new UsageConfiguration(prop);
  }

  public void test() throws Exception {
    Diagnostic.setLogStream();
    try (final TestDirectory dir = new TestDirectory()) {
      final File cfgFile = UsageConfiguration.createSimpleConfigurationFile(new File(dir, "filename"), dir.getPath(), "localhost:8080");
      final UsageConfiguration conf = new UsageConfiguration(cfgFile);
      assertEquals(dir.getPath(), conf.getUsageDir());
      assertEquals("localhost:8080", conf.getUsageHost());
    }
  }
}
