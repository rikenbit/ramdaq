//
//  Utility functions used in nf-core DSL2 module files
//

//
// Function to initialise default values and to generate a Groovy Map of available options for nf-core modules
//
def initOptions(Map args) {
    def Map options = [:]
    options.args            = args.args ?: ''
    options.args2           = args.args2 ?: ''
    options.args3           = args.args3 ?: ''
    options.publish_by_meta = args.publish_by_meta ?: []
    options.publish_dir     = args.publish_dir ?: ''
    options.publish_files   = args.publish_files
    options.suffix          = args.suffix ?: ''
    return options
}
