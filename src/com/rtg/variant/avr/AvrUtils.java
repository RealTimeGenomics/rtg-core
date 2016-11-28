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
package com.rtg.variant.avr;

import java.io.File;
import java.util.Arrays;

import com.rtg.util.Environment;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;

/**
 * Utility methods for locating AVR models.
 */
public final class AvrUtils {

  private AvrUtils() { }

  /** Property name that says where to load AVR models from */
  public static final String ENVIRONMENT_MODELS_DIR = "models.dir";
  /** flag for loading an AVR model to use for predictions */
  public static final String AVR_MODEL_FILE_FLAG = "avr-model";

  // Default AVR model name looked for within the above dir
  private static final String MODEL_DEFAULT = "illumina-exome.avr";
  // Name to explicitly specify no AVR model to be used
  private static final String AVR_NONE_NAME = "none";

  /**
   * Initialise flag for AVR model to use for prediction
   *
   * @param flags shared flags
   * @param anonymous if true, register the flag as an anonymous flag
   * @return the newly registered flag
   */
  public static Flag initAvrModel(final CFlags flags, boolean anonymous) {
    return initAvrModel(flags, anonymous, MODEL_DEFAULT);
  }

  /**
   * Initialise flag for AVR model to use for prediction
   *
   * @param flags shared flags
   * @param anonymous if true, register the flag as an anonymous flag
   * @param defaultModel specify the name of the model to be used as a default value
   * @return the newly registered flag
   */
  public static Flag initAvrModel(final CFlags flags, boolean anonymous, String defaultModel) {
    String description = "name of AVR model to use when scoring variants";
    File defaultModelFile = null;
    final String modelDirName = Environment.getEnvironmentMap().get(ENVIRONMENT_MODELS_DIR);
    if (modelDirName != null) {
      final File modelDirFile = new File(modelDirName);
      if (modelDirFile.exists() && modelDirFile.isDirectory()) {
        final String[] models = modelDirFile.list((dir, name) -> name.endsWith(".avr"));
        if (!anonymous && models != null && models.length > 0) {
          final String[] newModels = new String[models.length + 1];
          newModels[0] = AVR_NONE_NAME;
          System.arraycopy(models, 0, newModels, 1, models.length);
          Arrays.sort(newModels);
          description += ". Allowed values are " + Arrays.toString(newModels) + " or a path to a model file";
          if (Arrays.asList(newModels).contains(defaultModel)) {
            defaultModelFile = new File(defaultModel);
          }
        }
      }
    }
    final Flag modelFlag = anonymous
        ? flags.registerRequired(File.class, "MODEL", description).setCategory(CommonFlagCategories.REPORTING)
        : flags.registerOptional(AVR_MODEL_FILE_FLAG, File.class, "MODEL", description).setCategory(CommonFlagCategories.REPORTING);
    if (!anonymous && defaultModelFile != null) {
      modelFlag.setParameterDefault(defaultModelFile);
    }
    return modelFlag;
  }

  /**
   * Gets the AVR model to use for prediction. If the user has supplied a model name, it will be first
   * searched for as an actual file name, or secondly as a name within the environmental model directory
   * (if configured)
   *
   * @param flags shared flags
   * @param anonymous true if the model was registered as the only anonymous flag
   * @return a File pointing at the specified AVR model, or null if no AVR model is to be used
   * @throws InvalidParamsException if the user specified a model that does not exist, or a models directory which does
   * not exist or is not a directory.
   */
  public static File getAvrModel(final CFlags flags, boolean anonymous) {
    final String modelDirName = Environment.getEnvironmentMap().get(ENVIRONMENT_MODELS_DIR);
    File modelsDirFile = null;
    if (modelDirName != null) {
      modelsDirFile = new File(modelDirName);
      if (!modelsDirFile.exists() || !modelsDirFile.isDirectory()) {
        throw new InvalidParamsException("The AVR models directory cannot be found or is not a directory: " + modelDirName);
      }
    }

    final Flag flag = anonymous ? flags.getAnonymousFlag(0) : flags.getFlag(AVR_MODEL_FILE_FLAG);
    if (flag == null) {
      return null;
    }
    if (flag.getValue() != null) {
      File avrModel = (File) flag.getValue();
      if (AVR_NONE_NAME.equals(avrModel.toString())) {
        return null;
      } else if (avrModel.exists()) {
        return avrModel;
      } else {
        if (modelsDirFile != null) {
          avrModel = new File(modelsDirFile, avrModel.getName());
        }
        if (!avrModel.exists()) {
          throw new InvalidParamsException("The specified AVR model could not be found: " + flag.getValue());
        }
        return avrModel;
      }
    }
    return null;
  }
}
