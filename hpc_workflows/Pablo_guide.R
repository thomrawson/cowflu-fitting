## To run once:
# install.packages(
#   "hipercow",
#   repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))
library(hipercow)

## This will make a "hipercow" folder in your working directory that saves all relevant stuff
hipercow_init()

## This will configure the driver for Windows, there's definitely something for mac though, check with WEs!
hipercow_configure(driver = "windows")

## Check if it worked:
hipercow_configuration()

## Test task:
id <- task_create_expr(sessionInfo())
## See the submitted task's status:
task_status(id)
## See it's log:
task_log_show(id)
## Save any outputs:
output <- task_result(id)

################################################################
## When submitting tasks, you will need to log-in:
windows_authenticate()

## Loading required packages via conan can be done in a variety of ways,
## The easiest is to have a pkgdepends.txt file in your working directory
## See example on the hipercow website. Then, when it contains a list of your packages
## run:
hipercow_provision()

## You can list how many cores you want to ask for with this:
resources <- hipercow_resources(cores = 32)

##Then, within orderly, you just go ahead and submit a task like so:

task <- task_create_expr(orderly2::orderly_run("example_task"),
                         resources = resources)

task_status(task)
task_log_show(task)
task_info(task)$times
task_outputs <- task_result(task)

## If it's NOT orderly, and you instead have a custom function to run,
## Make a separate .R file that contains your function `my_function()` and then give
## hipercow this source of functions:
hipercow_environment_create(sources = "my_source_functions.R")

##Then you can submit those functions like above:
task2 <- task_create_expr(my_function(), resources = resources)
