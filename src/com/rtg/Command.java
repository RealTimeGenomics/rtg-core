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
package com.rtg;

import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Locale;

import com.rtg.launcher.AbstractCli;
import com.rtg.util.License;

/**
 * Wraps up a command that can be executed by name from a product entry and has license control, release status etc.
 */
public class Command {

  private final String mCommandName;
  private final CommandCategory mCategory;
  private final ReleaseLevel mReleaseLevel;
  private final String mLicenceProperty;
  private final String mAltLicenceProperty; // Used for backward compatibility with existing licenses when we change module names. May be null
  private final AbstractCli mCli;

  Command(final AbstractCli cli, final CommandCategory category, final ReleaseLevel level) {
    this(cli, category, level, null);
  }

  Command(final AbstractCli cli, final CommandCategory category, final ReleaseLevel level, final String altKey) {
    this(cli, cli.moduleName(), category, level, License.LICENSE_KEY_PREFIX + cli.moduleName(), altKey);
  }

  Command(final AbstractCli cli, final String commandName, final CommandCategory category, final ReleaseLevel level, final String key, final String altKey) {
    assert level != null;
    assert key != null;
    mCli = cli; // Should only be null if not implemented using the AbstractCli framework, and must therefore override mainInit
    mCommandName = commandName.toUpperCase(Locale.ROOT);
    mCategory = category;
    mReleaseLevel = level;
    mLicenceProperty = key.toLowerCase(Locale.ROOT);
    mAltLicenceProperty = altKey;
  }

  /**
   * Main init for running this module
   * @param args args for the module
   * @param out output stream to write to
   * @param err print stream to write to in case of error
   * @return integer return code
   */
  public int mainInit(final String[] args, final OutputStream out, final PrintStream err) {
    if (mCli == null) {
      throw new RuntimeException("Incorrectly configured module");
    }
    return mCli.mainInit(args, out, err);
  }

  /**
   * @return name of module
   */
  public String getCommandName() {
    return mCommandName;
  }

  @Override
  public String toString() {
    return mCommandName;
  }

  /**
   * @return key value used in license to determine if module is licensed
   */
  public String getLicenceKeyName() {
    return mLicenceProperty;
  }

  /**
   * @return an optional alternative license key name that the module may also be enabled by (for backwards compatibility).
   */
  public String getAltLicenceKeyName() {
    return mAltLicenceProperty;
  }

  /**
   * @return category of module
   */
  public CommandCategory getCategory() {
    return mCategory;
  }

  /**
   * @return release level of module
   */
  public ReleaseLevel getReleaseLevel() {
    return mReleaseLevel;
  }

  /**
   * @return true if should not be visible in help
   */
  public boolean isHidden() {
    return !(mReleaseLevel == ReleaseLevel.BETA || mReleaseLevel == ReleaseLevel.GA);
  }

  /**
   * @return true if current license may run this module
   */
  public boolean isLicensed() {
    return License.isPropertyLicensed(mLicenceProperty)
        || (mAltLicenceProperty != null && License.isPropertyLicensed(mAltLicenceProperty));
  }

  /**
   * @return string representation of license status
   */
  public String licenseStatus() {
    return mAltLicenceProperty != null && License.isPropertyLicensed(mAltLicenceProperty)
        ? License.getExpiryStatus(mAltLicenceProperty)
        : License.getExpiryStatus(mLicenceProperty);
  }

}
