CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

# Setting policies, to avoid warnings
CMAKE_POLICY( VERSION 2.8 )
IF( POLICY CMP0028 )
	CMAKE_POLICY( SET CMP0028 OLD )
ENDIF( POLICY CMP0028 )
IF( POLICY CMP0043 )
	CMAKE_POLICY( SET CMP0043 OLD )
ENDIF( POLICY CMP0043 )
IF( POLICY CMP0046 )
	CMAKE_POLICY( SET CMP0046 OLD )
ENDIF( POLICY CMP0046 )

PROJECT( CartPoleExperiment )

### Enabling testing
INCLUDE( CTest )
ENABLE_TESTING()

### Using FARSA for generating the plugin
FIND_PACKAGE( FARSA REQUIRED )
SET( FARSA_USE_EXPERIMENTS TRUE )
INCLUDE( ${FARSA_USE_FILE} )

# We have to increase the minimum required version if using QT5
IF( FARSA_USE_QT5 )
	CMAKE_MINIMUM_REQUIRED( VERSION 2.8.9 )

	# Setting policies, to avoid warnings
	CMAKE_POLICY( VERSION 2.8.9 )
	IF( POLICY CMP0028 )
		CMAKE_POLICY( SET CMP0028 OLD )
	ENDIF( POLICY CMP0028 )
	IF( POLICY CMP0043 )
		CMAKE_POLICY( SET CMP0043 OLD )
	ENDIF( POLICY CMP0043 )
	IF( POLICY CMP0046 )
		CMAKE_POLICY( SET CMP0046 OLD )
	ENDIF( POLICY CMP0046 )
ENDIF( FARSA_USE_QT5 )

### set sources and headers files
FILE( GLOB SW_SRCS ./src/*.cpp )
FILE( GLOB SW_HDRS ./include/*.h )
### Add headers to sources for moc'ing up as well as resources
IF(FARSA_USE_QT5)
	#QT5_ADD_RESOURCES(SW_SRCS "./experimentTextures.qrc")
ELSE(FARSA_USE_QT5)
	QT4_WRAP_CPP(SW_SRCS ${SW_HDRS} OPTIONS "-nw")
	#QT4_ADD_RESOURCES(SW_SRCS "./experimentTextures.qrc")
ENDIF(FARSA_USE_QT5)

INCLUDE_DIRECTORIES( ./include )

### Use the ADD_FARSAPLUGIN for creating the target
### The plugin is a shader library loaded at runtime by total99
ADD_FARSAPLUGIN( ${PROJECT_NAME} ${SW_SRCS} ${SW_HDRS} )

# Also installing the plugin configuration
FARSA_PLUGIN_INSTALL_CONF( ${PROJECT_NAME} "conf" )

### Test-related stuffs
IF( BUILD_TESTING )
	### Setting the build name to a default value in case of tests
	SET( BUILDNAME "${PROJECT_NAME}" CACHE STRING "Name of build on the dashboard" )
	MARK_AS_ADVANCED( BUILDNAME )

	### Copying the file to configure CTest to the build directory
	CONFIGURE_FILE( CTestCustom.cmake ${CMAKE_BINARY_DIR}/CTestCustom.cmake COPYONLY )

	### Adding the test to install the plugin
	ADD_FARSAPLUGIN_INSTALL_TEST()
ENDIF( BUILD_TESTING )

