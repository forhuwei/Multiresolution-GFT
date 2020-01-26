function QStep = QStep_Compute(QP)

QStep_base = [0.625 0.6875 0.8125 0.875 1 1.125 1.25];
% QStep_base = [0.625 0.702 0.787 0.884 0.992 1.114 1.25];

if QP <= 6
    QStep = QStep_base(QP+1);
    return
end

QStep = QStep_base(mod(QP,6)+1)*2^floor(QP/6);