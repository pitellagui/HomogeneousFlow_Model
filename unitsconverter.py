"""
Módulo para conversão de unidades de temperatura e pressão.

Métodos:
----------
    Temperature(value, from_unity, to_unity): Converte valores de temperatura
    Pressure(value, from_unity, to_unity): Converte valores de pressão

Parâmetros:
----------
        value (float): Valor numérico da temperatura a ser convertida
        from_unity (str): Unidade de origem ('C', 'F', 'K' ou 'R')
        to_unity (str): Unidade de destino ('C', 'F', 'K' ou 'R')

Retorna:
--------
        float: Valor convertido na unidade desejada

    Escalas suportadas:
    - Temperatura:
        - Celsius (°C)
        - Fahrenheit (°F)
        - Kelvin (K)
        - Rankine (°R)
    - Pressão:
        - psi (libra-força por polegada quadrada)
        - psia (psi absoluto - mesmo valor numérico que psi)
        - Pa (Pascal)
        - bar
        - atm (atmosfera padrão)
    """


def Temperature(value, from_unity, to_unity):

    if from_unity == 'C':
        if to_unity == 'F':
            return (value * 9 / 5) + 32
        elif to_unity == 'K':
            return value + 273.15
        elif to_unity == 'R':
            return (value + 273.15) * 9 / 5


    elif from_unity == 'F':
        if to_unity == 'C':
            return (value - 32) * 5 / 9
        elif to_unity == 'K':
            return (value - 32) * 5 / 9 + 273.15
        elif to_unity == 'R':
            return value + 459.67


    elif from_unity == 'K':
        if to_unity == 'C':
            return value - 273.15
        elif to_unity == 'F':
            return (value - 273.15) * 9 / 5 + 32
        elif to_unity == 'R':
            return value * 9 / 5


    elif from_unity == 'R':
        if to_unity == 'C':
            return (value - 491.67) * 5 / 9
        elif to_unity == 'F':
            return value - 459.67
        elif to_unity == 'K':
            return value * 5 / 9

    if from_unity == to_unity:
        return value

    raise ValueError(f"Conversão de {from_unity} para {to_unity} não suportada")


def Pressure(value, from_unity, to_unity):

    if to_unity == 'psi':
        if from_unity == 'Pa':
            return value / 6894.76
        elif from_unity == 'bar':
            return value * 14.5038
        elif from_unity == 'atm':
            return value * 14.6959
        elif from_unity == 'psi':
            return value


    elif from_unity == 'psi':
        if to_unity == 'Pa':
            return value * 6894.76
        elif to_unity == 'bar':
            return value / 14.5038
        elif to_unity == 'atm':
            return value / 14.6959


    elif from_unity == 'bar' and to_unity == 'atm':
        return value * 0.986923

    if from_unity == to_unity:
        return value

    raise ValueError(f"Conversão de {from_unity} para {to_unity} não suportada")