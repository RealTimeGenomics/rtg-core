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
import java.util.ArrayList;
import java.util.Locale;
import java.util.Properties;
import java.util.UUID;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Delegates usage tracking messages to appropriate recording implementation. All external usage tracking
 * calls should go through here.
 */
public final class UsageTracking {

  static final String SEE_MANUAL = "Please consult the user manual section \"Advanced Installation Configuration\".";

  enum UsageDestination {
    NONE,
    FILE_OR_SERVER,
    SERVER_ONLY
  }

  // Name of usage tracking properties obtained from license
  static final String REQUIRE_USAGE = "require_usage";
  static final String USAGE_DESTINATION = "usage_destination";

  /** Records messages generated to the usage clients so can test if the calls are being done. */
  private final ArrayList<String> mUsageMessages = new ArrayList<>();

  private final boolean mRequireUsage;
  private final UsageDestination mUsageDestination;
  private final String mModuleName;
  private final UUID mRunId;
  private final UsageTrackingClient mClient;
  private final UsageConfiguration mUsageConfiguration;


  /**
   * Sets up usage tracking using default configuration path
   * @param properties contains license configured properties (whether usage is required, and which usage tracking modes are allowed)
   * @param modulename name of module being run
   * @param runId unique ID for the run
   * @param suppress true if we should ignore all usage messages
   * @throws IOException when an IO error occurs
   */
  public UsageTracking(Properties properties, String modulename, UUID runId, boolean suppress) throws IOException {
    this(properties, modulename, runId, null, suppress);
  }

  /**
   * Sets up usage tracking using a specific configuration path
   * @param licenseProperties contains license configured properties (whether usage is required, and which usage tracking modes are allowed)
   * @param modulename name of module being run
   * @param runId unique ID for the run
   * @param configFileOverride path of file to load install specific usage tracking configuration from (directory for file based tracking, or host for server based tracking)
   * @param suppress true if we should use the null reporter (effectively ignoring all usage messages)
   * @throws IOException if an IO error occurs
   */
  UsageTracking(Properties licenseProperties, String modulename, UUID runId, File configFileOverride, boolean suppress) throws IOException {
    final String usageValue = licenseProperties.getProperty(REQUIRE_USAGE);
    mRequireUsage = usageValue != null && Boolean.valueOf(usageValue);
    final String destinationValue = licenseProperties.getProperty(USAGE_DESTINATION);
    mUsageDestination = destinationValue == null ? UsageDestination.NONE : UsageDestination.valueOf(destinationValue.toUpperCase(Locale.ROOT));
    mUsageConfiguration = configFileOverride == null ? new UsageConfiguration() : new UsageConfiguration(configFileOverride);
    mModuleName = modulename;
    mRunId = runId;
    if (suppress) {
      mClient = new NullUsageTrackingClient();
    } else {
      if (mUsageConfiguration.isEnabled() && mUsageConfiguration.getUsageHost() != null) {
        mClient = new HttpUsageTrackingClient(mUsageConfiguration.getUsageHost(), mUsageConfiguration, requireUsage());
      } else if (mUsageConfiguration.isEnabled() && allowFileTracking() && mUsageConfiguration.getUsageDir() != null) {
        mClient = new FileUsageTrackingClient(new File(mUsageConfiguration.getUsageDir()), mUsageConfiguration, requireUsage());
      } else if (requireUsage()) {
        throw new NoTalkbackSlimException("Usage logging is required by license, but has not been correctly configured during install. " + SEE_MANUAL);
      } else {
        mClient = new NullUsageTrackingClient();
      }
    }
  }

  /**
   * @return the usage configuration
   */
  public UsageConfiguration getUsageConfiguration() {
    return mUsageConfiguration;
  }

  /**
   * @return a String with a human readable version of all the usage messages so far.
   */
  public String usageLog() {
    return mUsageMessages.toString();
  }

  /**
   * records a start usage tracking message to whichever tracking endpoint is configured
   */
  public void recordBeginning() {
    final String msg = "Usage beginning module=" + mModuleName + " runId=" + mRunId;
    mUsageMessages.add(msg);
    Diagnostic.developerLog(msg);
    mClient.recordBeginning(mModuleName, mRunId);
  }

  /**
   * records an end usage tracking message to whichever tracking endpoint is configured
   * @param metric module specific usage metric
   * @param success true if run succeeded
   */
  public void recordEnd(long metric, boolean success) {
    final String msg = "Usage end module=" + mModuleName + " runId=" + mRunId + " metric=" + metric + " success=" + success;
    mUsageMessages.add(msg);
    Diagnostic.developerLog(msg);
    mClient.recordEnd(metric, mModuleName, mRunId, success);
  }

  /**
   * if true the {@link NullUsageTrackingClient} is not allowed
   * @return if a usage tracking client must be configured
   */
  boolean requireUsage() {
    return mRequireUsage;
  }

  /**
   * if true the {@link FileUsageTrackingClient} is allowed
   * @return if file based usage tracking is allowed
   */
  boolean allowFileTracking() {
    return mUsageDestination == UsageDestination.FILE_OR_SERVER;
  }
}
