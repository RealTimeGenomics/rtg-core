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
package com.rtg.util;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

import com.rtg.util.diagnostic.ErrorType;


/**
 * Utility class for property file loading and parsing.
 */
public final class PropertiesUtils {
  private PropertiesUtils() { }

  static final String VARIANT_PRIORS = "com/rtg/variant/priors/";
  static final String ALIGNMENT_PENALTIES = "com/rtg/alignment/penalties/";

  /**
   * Types of properties files, and where they live
   */
  public static enum PropertyType {
    /** Error property sub-directory **/
   ERROR_PROPERTY(VARIANT_PRIORS + "error/"),
    /** Prior property sub-directory **/
    PRIOR_PROPERTY(VARIANT_PRIORS + "prior/"),
    /** CNV property sub-directory **/
    CNV_PROPERTY(VARIANT_PRIORS + "cnv/"),
    /** Alignment error penalties */
    ALIGNMENT_PROPERTY_TYPE(ALIGNMENT_PENALTIES);


    final String mPath;
    private PropertyType(String path) {
      mPath = path;
    }

    /**
     * @return the path to the properties
     */
    public String path() {
      return mPath;
    }
  }

  /**
   * Method for loading a priors properties resource file or a properties file.
   * @param propsFile the priors property file or resource to load.
   * @param propertyType the type of property file being loaded (sub-directory).
   * @return Loaded properties object.
   * @throws InvalidParamsException when the property file does not exist.
   * @throws IOException when there is an error reading the property files.
   */
  public static Properties getPriorsResource(final String propsFile, final PropertyType propertyType) throws InvalidParamsException, IOException {
    final Properties pr = new Properties();
    InputStream props = null;
    try {
      final InputStream def = Resources.getResourceAsStream(propertyType.path() + propsFile + ".properties");
      if (def != null) {
        props = def;
      } else {
        try {
          props = new FileInputStream(propsFile);
        } catch (final FileNotFoundException e) {
          throw new InvalidParamsException(ErrorType.INFO_ERROR, "Invalid prior option \"" + propsFile + "\"");
        }
      }
      try {
        pr.load(props);
      } catch (final IOException e) {
        throw new InvalidParamsException(ErrorType.PROPS_LOAD_FAILED, "Priors", propsFile, e.getMessage());
      } catch (final IllegalArgumentException e) {
        throw new InvalidParamsException(ErrorType.PROPS_INVALID, "Priors", propsFile);
      }
    } finally {
      if (props != null) {
        props.close();
      }
    }
    return pr;
  }
}
