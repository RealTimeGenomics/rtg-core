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
