/*
 * Copyright (c) 2018. Real Time Genomics Limited.
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

import com.rtg.launcher.AbstractCli;
import com.rtg.util.License;

/**
 * Command for utilities that are too simple to bother with license restrictions
 */
public class LicensedCommand extends Command {

  LicensedCommand(AbstractCli cli, final CommandCategory category, final ReleaseLevel level) {
    super(cli, cli.moduleName(), cli.description(), category, level, License.LICENSE_KEY_PREFIX + cli.moduleName(), null);
  }

  @Override
  public boolean isLicensed() {
    return License.isPropertyLicensed(getLicenceKeyName())
      || (getAltLicenceKeyName() != null && License.isPropertyLicensed(getAltLicenceKeyName()));
  }

  @Override
  public String licenseStatus() {
    return getAltLicenceKeyName() != null && License.isPropertyLicensed(getAltLicenceKeyName())
      ? License.getExpiryStatus(getAltLicenceKeyName())
      : License.getExpiryStatus(getLicenceKeyName());
  }

}
