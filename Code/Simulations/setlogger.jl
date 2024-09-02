function setlogger(logfile)
    FileLogger(logfile, append=true) |> global_logger
end


# function setlogger()
#     TeeLogger(
#         MinLevelLogger(FileLogger("/tmpinfo.log"), Logging.Info),
#         MinLevelLogger(FileLogger("/tmp/warn.log"), Logging.Warn),
#     ) |> global_logger
# end 



# function setlogger()
#     # write each process' logs to its own file
#     TeeLogger(
#         MinLevelLogger(FileLogger("/tmpinfo_$(my_id()).log"), Logging.Info),
#         MinLevelLogger(FileLogger("/tmp/warn_$(my_id()).log"), Logging.Warn),
#     ) |> global_logger
# end 