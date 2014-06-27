Config::Config
{
    // Read entire configuration file as is
    readConfigFile(configFileName);

    // Update user values with command-line values
    setUserSettings(commandLineParameters);
    
    const std::string& configType = userValue(CONFIG_TYPE_NAME);

    if (configType == "speciationextinction") {
        initializeDiversification();
    } else if (configType == "trait") {
        initializeTrait();
    } else {
        // throw an error
    }

    validateUserSettings();
}
