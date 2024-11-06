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

import static com.rtg.variant.avr.AbstractModelBuilder.MODEL_AVR_VERSION;
import static com.rtg.variant.avr.AbstractModelBuilder.MODEL_PROPERTY_COMMAND_LINE;
import static com.rtg.variant.avr.AbstractModelBuilder.MODEL_PROPERTY_DATE;
import static com.rtg.variant.avr.AbstractModelBuilder.MODEL_PROPERTY_DERIVED_ANNOTATIONS;
import static com.rtg.variant.avr.AbstractModelBuilder.MODEL_PROPERTY_FORMAT_ANNOTATIONS;
import static com.rtg.variant.avr.AbstractModelBuilder.MODEL_PROPERTY_INFO_ANNOTATIONS;
import static com.rtg.variant.avr.AbstractModelBuilder.MODEL_PROPERTY_MODEL_ID;
import static com.rtg.variant.avr.AbstractModelBuilder.MODEL_PROPERTY_QUAL_ANNOTATION;
import static com.rtg.variant.avr.AbstractModelBuilder.MODEL_PROPERTY_TYPE;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Properties;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.util.Pair;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.LineWriter;

/**
 */
public class AvrStatsCli extends AbstractCli {

  private static final String DUMP_MODEL_FLAG = "Xdump-model";
  private static final String DUMP_PROPERTIES_FLAG = "Xdump-properties";
  private static final String UPGRADE_MODEL_FLAG = "Xupgrade-model";

  @Override
  public String moduleName() {
    return "avrstats";
  }

  @Override
  public String description() {
    return "print statistics about an AVR model";
  }

  @Override
  protected void initFlags() {
    CommonFlagCategories.setCategories(mFlags);
    mFlags.setDescription("Print statistics that describe an AVR model.");
    mFlags.registerOptional(DUMP_PROPERTIES_FLAG, "if set, output the raw model properties").setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerOptional(DUMP_MODEL_FLAG, "if set, output a verbose representation of the model").setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerOptional(UPGRADE_MODEL_FLAG, File.class, CommonFlags.FILE, "if set, re-save the model (to upgrade model version)").setCategory(CommonFlagCategories.UTILITY);

    AvrUtils.initAvrModel(mFlags, true);
  }

  private final ArrayList<Pair<String, String>> mProperties = new ArrayList<>();
  {
    mProperties.add(new Pair<>("Date built", MODEL_PROPERTY_DATE));
    mProperties.add(new Pair<>("AVR Version", MODEL_AVR_VERSION));
    mProperties.add(new Pair<>("AVR-ID", MODEL_PROPERTY_MODEL_ID));
    mProperties.add(new Pair<>("Type", MODEL_PROPERTY_TYPE));
    mProperties.add(new Pair<>("QUAL used", MODEL_PROPERTY_QUAL_ANNOTATION));
    mProperties.add(new Pair<>("INFO fields", MODEL_PROPERTY_INFO_ANNOTATIONS));
    mProperties.add(new Pair<>("FORMAT fields", MODEL_PROPERTY_FORMAT_ANNOTATIONS));
    mProperties.add(new Pair<>("Derived fields", MODEL_PROPERTY_DERIVED_ANNOTATIONS));
    mProperties.add(new Pair<>("Parameters", MODEL_PROPERTY_COMMAND_LINE));
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final File modelFile;
    modelFile = AvrUtils.getAvrModel(mFlags, true);
    if (modelFile == null) {
      throw new NoTalkbackSlimException("No model file specified and no default model available.");
    }
    try (LineWriter lw = new LineWriter(new OutputStreamWriter(out))) {
      final ModelFactory fact = new ModelFactory(modelFile);
      final Properties props = fact.getModelProperties();
      int maxLen = "Location".length();
      for (Pair<String, String> field : mProperties) {
        maxLen = Math.max(maxLen, field.getA().length());
      }
      lw.writeln("Location" + StringUtils.spaces(maxLen - 8) + ": " + modelFile);
      for (Pair<String, String> field : mProperties) {
        final String label = field.getA();
        final String value = props.getProperty(field.getB());
        if (value != null && value.length() > 0) {
          lw.writeln(field.getA() + StringUtils.spaces(maxLen - label.length()) + ": " + value);
        }
      }
      lw.writeln();

      if (mFlags.isSet(DUMP_PROPERTIES_FLAG)) {
        lw.writeln("Properties:");
        for (String property : props.stringPropertyNames()) {
          lw.writeln(property + "=" + props.getProperty(property));
        }
        lw.writeln();
      }

      if (mFlags.isSet(DUMP_MODEL_FLAG)) {
        lw.writeln("Full Model:");
        lw.writeln(fact.getModel().toString());
        lw.writeln();
      }

      if (mFlags.isSet(UPGRADE_MODEL_FLAG)) {
        final File outFile = (File) mFlags.getValue(UPGRADE_MODEL_FLAG);
        lw.writeln("Writing Model To: " + outFile);
        AbstractModelBuilder.save(outFile, fact.getModelProperties(), fact.getModel());
      }
    }
    return 0;
  }
}
