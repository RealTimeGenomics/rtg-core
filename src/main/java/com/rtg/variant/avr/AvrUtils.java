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
  private static final String MODEL_DEFAULT = "illumina-wgs.avr";
  // Name to explicitly specify no AVR model to be used
  private static final String AVR_NONE_NAME = "none";

  /**
   * Initialise flag for AVR model to use for prediction
   *
   * @param flags shared flags
   * @param anonymous if true, register the flag as an anonymous flag
   * @return the newly registered flag
   */
  public static Flag<File> initAvrModel(final CFlags flags, boolean anonymous) {
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
  public static Flag<File> initAvrModel(final CFlags flags, boolean anonymous, String defaultModel) {
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
      } else {
        throw new InvalidParamsException("Invalid installation - the AVR models directory cannot be found or is not a directory: " + modelDirName);
      }
    }
    final Flag<File> modelFlag = anonymous
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
        throw new InvalidParamsException("Invalid installation - the AVR models directory cannot be found or is not a directory: " + modelDirName);
      }
    }

    final Flag<?> flag = anonymous ? flags.getAnonymousFlag(0) : flags.getFlag(AVR_MODEL_FILE_FLAG);
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
          avrModel = new File(modelsDirFile, avrModel.getPath());
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
