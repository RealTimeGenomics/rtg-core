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
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Map;
import java.util.Properties;

import com.rtg.util.Environment;
import com.rtg.util.io.FileUtils;

/**
 */
public class UsageConfiguration {


  static final String DEFAULT_USAGE_HOST = "https://api.realtimegenomics.com/usage/submit.php";

  static final String ENABLE_USAGE = "usage";
  static final String USAGE_DIR = "usage.dir";
  static final String USAGE_HOST = "usage.host";
  static final String USAGE_LOG_USERNAME = "usage.log.username";
  static final String USAGE_LOG_HOSTNAME = "usage.log.hostname";
  static final String USAGE_LOG_COMMANDLINE = "usage.log.commandline";

  private final Properties mProp;

  /**
   * Construct default usage configuration from environment-supplied settings
   */
  UsageConfiguration() {
    mProp = new Properties();
    final Map<String, String> env = Environment.getEnvironmentMap();
    final Boolean enabled = Boolean.valueOf(env.get(ENABLE_USAGE));
    mProp.setProperty(ENABLE_USAGE, enabled.toString());
    if (enabled) {
      // Pull the properties we care about out of system properties
      for (String property : new String[] {USAGE_DIR, USAGE_HOST, USAGE_LOG_USERNAME, USAGE_LOG_HOSTNAME, USAGE_LOG_COMMANDLINE}) {
        if (env.get(property) != null) {
          mProp.setProperty(property, env.get(property));
        }
      }

      // Set default usage host if other destinations have not been set
      if (env.get(USAGE_DIR) == null && env.get(USAGE_HOST) == null) {
        mProp.setProperty(USAGE_HOST, DEFAULT_USAGE_HOST);
      }
    }
  }

  /**
   * Load usage configuration from an external properties file
   * @param configFile contains configuration settings
   * @throws IOException if the configuration file could not be read
   */
  UsageConfiguration(File configFile) throws IOException {
    mProp = new Properties();
    try (InputStream is = FileUtils.createInputStream(configFile, false)) {
      mProp.load(is);
    }
    mProp.setProperty(ENABLE_USAGE, "true"); // If you're bothering to use an explicit configuration file, it's enabled
  }

  UsageConfiguration(Properties prop) {
    mProp = prop;
    mProp.setProperty(ENABLE_USAGE, "true"); // If you're bothering to use an explicit configuration file, it's enabled
  }

  static File createSimpleConfigurationFile(File propFile, String usageDir, String usageHost) throws IOException {
    final Properties prop = new Properties();
    if (usageDir != null) {
      prop.setProperty(USAGE_DIR, usageDir);
    }
    if (usageHost != null) {
      prop.setProperty(USAGE_HOST, usageHost);
    }
    try (OutputStream os = FileUtils.createOutputStream(propFile)) {
      prop.store(os, "");
    }
    return propFile;
  }

  public boolean isEnabled() {
    return Boolean.valueOf(mProp.getProperty(ENABLE_USAGE));
  }

  public String getUsageDir() {
    return mProp.getProperty(USAGE_DIR);
  }

  public String getUsageHost() {
    return mProp.getProperty(USAGE_HOST);
  }

  /**
   * @return true if usage messages should include the current user
   */
  public boolean logUsername() {
    return Boolean.valueOf(mProp.getProperty(USAGE_LOG_USERNAME, "false"));
  }

  /**
   * @return true if usage messages should include the current machine name
   */
  public boolean logHostname() {
    return Boolean.valueOf(mProp.getProperty(USAGE_LOG_HOSTNAME, "false"));
  }

  /**
   * @return true if usage messages should include the current command line
   */
  public boolean logCommandLine() {
    return Boolean.valueOf(mProp.getProperty(USAGE_LOG_COMMANDLINE, "false"));
  }
}
