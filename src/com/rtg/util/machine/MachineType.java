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

package com.rtg.util.machine;

import java.util.Locale;

import com.rtg.util.EnumHelper;
import com.rtg.util.PseudoEnum;

/**
 * Enumeration of different machine types and associated information.
 */
public final class MachineType implements PseudoEnum {

  private static int sSequenceNumber = -1;

  /** The platform type for Illumina, regardless of whether it is paired or single end */
  public static final String PLAT_ILLUMINA = "ILLUMINA";
  /** The platform type for 454, regardless of whether it is paired or single end */
  public static final String PLAT_454 = "LS454";

  /** Illumina single end. */
  public static final MachineType ILLUMINA_SE = new MachineType(++sSequenceNumber, "illumina_se", "illumina", PLAT_ILLUMINA, null);

  /** Illumina paired end. */
  public static final MachineType ILLUMINA_PE = new MachineType(++sSequenceNumber, "illumina_pe", "illumina", PLAT_ILLUMINA, MachineOrientation.FR);

  /** Complete Genomics (paired end). */
  public static final MachineType COMPLETE_GENOMICS = new MachineType(++sSequenceNumber, "complete_genomics", "complete", "COMPLETE", MachineOrientation.TANDEM);

  /** Four Five Four paired end. */
  public static final MachineType FOURFIVEFOUR_PE = new MachineType(++sSequenceNumber, "454_pe", "ls454_pe", PLAT_454, null);

  /** Four Five Four single end. */
  public static final MachineType FOURFIVEFOUR_SE = new MachineType(++sSequenceNumber, "454_se", "ls454_se", PLAT_454, null);

  /** Ion Torrent (single end). */
  public static final MachineType IONTORRENT = new MachineType(++sSequenceNumber, "iontorrent", "iontorrent", "IONTORRENT", null);

  static final EnumHelper<MachineType> HELPER = new EnumHelper<>(MachineType.class, new MachineType[] {ILLUMINA_SE, ILLUMINA_PE, COMPLETE_GENOMICS, FOURFIVEFOUR_PE, FOURFIVEFOUR_SE, IONTORRENT});

  /**
   * @return list of the enum names
   */
  public static String[] names() {
    return HELPER.names();
  }

  /**
   * see {@link java.lang.Enum#valueOf(Class, String)}
   * @param str the name of the enum
   * @return the enum value
   */
  public static MachineType valueOf(final String str) {
    return HELPER.valueOf(str.toLowerCase(Locale.ROOT));
  }


  private final String mName;
  private final int mOrdinal;
  private final String mPriors;
  private final MachineOrientation mOrientation;
  private final String mPlatform;

  private MachineType(int ordinal, String name, String defaultPriors, String platform, MachineOrientation orientation) {
    mName = name;
    mOrdinal = ordinal;
    mPriors = defaultPriors;
    mOrientation = orientation;
    mPlatform = platform;
  }

  @Override
  public String toString() {
    return mName;
  }

  @Override
  public String name() {
    return mName;
  }

  @Override
  public int ordinal() {
    return mOrdinal;
  }

  /**
   * Gets the name of the default priors for this machine
   * @return the default priors.
   */
  public String priors() {
    return mPriors;
  }

  /**
   * Get the <code>MachineOrientation</code>.
   * @return the <code>MachineOrientation</code> (may be null if unknown or single end data).
   */
  public MachineOrientation orientation() {
    return mOrientation;
  }

  /**
   * @return the string defining the platform of this machine type
   */
  public String platform() {
    return mPlatform;
  }

  /**
   * Compare if the platforms are compatible, e.g. <code>ILLUMINA</code> is equal to <code>illumina</code> or <code>Illumina</code>
   * @param platform other platform
   * @return if the two platforms are compatible
   */
  public boolean compatiblePlatform(final String platform) {
    return mPlatform.equalsIgnoreCase(platform);
  }

  /**
   * {@link EnumHelper#values()}
   * @return as in link
   */
  public static MachineType[] values() {
    return HELPER.values();
  }
}
