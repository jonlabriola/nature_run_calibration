#python distance_correlation.py 25 uinterp uinterp 44  &
#python distance_correlation.py 25 uinterp uinterp 55  &
#python distance_correlation.py 25 uinterp uinterp 64  &
#python distance_correlation.py 25 uinterp uinterp 25  &

#python distance_correlation_pert.py 25 uinterp uinterp 44  &
#python distance_correlation_pert.py 25 uinterp uinterp 55  &
#python distance_correlation_pert.py 25 uinterp uinterp 64  &
#python distance_correlation_pert.py 25 uinterp uinterp 25  &
#wait
python distance_correlation.py 25 uinterp uinterp 44 --vertical_localization  &
python distance_correlation.py 25 uinterp uinterp 55 --vertical_localization &
python distance_correlation.py 25 uinterp uinterp 64 --vertical_localization &
python distance_correlation.py 25 uinterp uinterp 25 --vertical_localization  &

python distance_correlation_pert.py 25 uinterp uinterp 44  --vertical_localization &
python distance_correlation_pert.py 25 uinterp uinterp 55  --vertical_localization &
python distance_correlation_pert.py 25 uinterp uinterp 64  --vertical_localization &
python distance_correlation_pert.py 25 uinterp uinterp 25  --vertical_localization &
wait
python distance_correlation.py 25 vinterp vinterp 44 --vertical_localization  &
python distance_correlation.py 25 vinterp vinterp 55 --vertical_localization &
python distance_correlation.py 25 vinterp vinterp 64 --vertical_localization &
python distance_correlation.py 25 vinterp vinterp 25 --vertical_localization  &

python distance_correlation_pert.py 25 vinterp vinterp 44  --vertical_localization &
python distance_correlation_pert.py 25 vinterp vinterp 55  --vertical_localization &
python distance_correlation_pert.py 25 vinterp vinterp 64  --vertical_localization &
python distance_correlation_pert.py 25 vinterp vinterp 25  --vertical_localization &
wait
python distance_correlation.py 25 qv qv 44 --vertical_localization  &
python distance_correlation.py 25 qv qv 55 --vertical_localization &
python distance_correlation.py 25 qv qv 64 --vertical_localization &
python distance_correlation.py 25 qv qv 25 --vertical_localization  &

python distance_correlation_pert.py 25 qv qv 44  --vertical_localization &
python distance_correlation_pert.py 25 qv qv 55  --vertical_localization &
python distance_correlation_pert.py 25 qv qv 64  --vertical_localization &
python distance_correlation_pert.py 25 qv qv 25  --vertical_localization &
wait
python distance_correlation.py 25 uinterp winterp 44 --vertical_localization  &
python distance_correlation.py 25 uinterp winterp 55 --vertical_localization &
python distance_correlation.py 25 uinterp winterp 64 --vertical_localization &
python distance_correlation.py 25 uinterp winterp 25 --vertical_localization  &

python distance_correlation_pert.py 25 uinterp winterp 44  --vertical_localization &
python distance_correlation_pert.py 25 uinterp winterp 55  --vertical_localization &
python distance_correlation_pert.py 25 uinterp winterp 64  --vertical_localization &
python distance_correlation_pert.py 25 uinterp winterp 25  --vertical_localization &
wait

