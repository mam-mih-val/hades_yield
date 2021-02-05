add_executable(yield yield.cc yield.h)
target_link_libraries(yield PUBLIC ROOT::EG at_task_main)

add_executable(acceptance acceptance.cc)
target_link_libraries(acceptance PUBLIC ROOT::EG at_task_main)