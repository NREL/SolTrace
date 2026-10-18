# Load the pysoltrace api from the parent directory ---
import ctypes, os, sys
sys.path.insert(1, os.path.join(sys.path[0], '../..'))

from pysoltrace import api, dot_h, found_in

from pysoltrace.timer import timer

t = timer()

if __name__ == '__main__':
    stapi = api()
    # print(f'\n\n\n\n{found_in(dot_h)}\n\n\n\n')
    # print(getattr(dot_h, 'args_st_write_results_csv'))
    # test = dot_h.args_st_write_results_csv()
    # print(test._fields_)
    # for f in test._fields_:
    #     print(getattr(test, f[0]))

    # test_payload = getattr(dot_h.st_api_call_args(), 'payload')
    # print(test_payload._fields_)
    # for f in test_payload._fields_:
    #     print(getattr(test_payload, f[0])._fields_)
    #     for pf in getattr(test_payload, f[0])._fields_:
    #         print(getattr(getattr(test_payload, f[0]), pf[0]))
    # dummy_array = (ctypes.c_double * 0)(*[])
    # test = stapi.generate_api_call(dot_h.st_api_call.CALL_ST_ADD_SUN,
    #                               ctypes.pointer(dot_h.args_sun(0, 608, 303, 1000, 5, 'g'.encode())),
    #                               ctypes.pointer(ctypes.c_double()),
    #                               ctypes.pointer(ctypes.c_double()))

    # # test = stapi.generate_api_call(dot_h.st_api_call.CALL_ST_SIM_SETUP, dot_h.st_runner_type_t.OPTIX)
    # print(test)
    # print(test.payload._fields_[test.type])
    # arg_name = test.payload._fields_[test.type][0]

    # for field in getattr(test.payload, arg_name)._fields_:
    #     temp_payload_attr = getattr(test.payload, test.payload._fields_[test.type][0])
    #     print(getattr(temp_payload_attr, field[0]))
    # # print(dummy_list)


    print('timing comparison')
    for _ in range(10):
        t.ic('NO batch')
        stapi.data.json.load('../sample.json')
        count = len(stapi.data.element)
        stapi.runner.setup(dot_h.st_runner_type_t.OPTIX)
        stapi.runner.run()
        stapi.runner.report()
        # stapi.write_results_csv('../batch_sample.csv')
        t.oc('NO batch')

    for _ in range(10):
        t.ic('batch')

        stapi.batch.data.json.load('../sample.json')
        rec = stapi.batch.data.element.num()
        stapi.batch.runner.setup(dot_h.st_runner_type_t.OPTIX)
        stapi.batch.runner.run()
        stapi.batch.runner.report()

        stapi.batch()
        pcount = ctypes.c_uint64()
        # stapi.batch([
        #     stapi.generate_api_call(dot_h.st_api_call.CALL_ST_READ_INPUT_JSON, f.read()),
        #     stapi.generate_api_call(dot_h.st_api_call.CALL_ST_NUM_ELEMENTS, ctypes.pointer(pcount)),
        #     stapi.generate_api_call(dot_h.st_api_call.CALL_ST_SIM_SETUP, dot_h.st_runner_type_t.OPTIX),
        #     stapi.generate_api_call(dot_h.st_api_call.CALL_ST_SIM_RUN_V2),
        #     stapi.generate_api_call(dot_h.st_api_call.CALL_ST_SIM_REPORT, 0),
        #     # stapi.generate_api_call(dot_h.st_api_call.CALL_ST_WRITE_RESULTS_CSV, b'../batch_sample.csv', 12)
        # ], True)
        t.oc('batch')

    print(count)
    print(rec.value)
    print(t)