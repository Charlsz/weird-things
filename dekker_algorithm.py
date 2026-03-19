import threading
import time

# Variables compartidas (Memoria global del sistema)
inventario = 100
flag = [False, False]  # Intención de entrar a la sección crítica
turn = 0               # Mecanismo de desempate

def modificar_inventario(mi_id):
    global inventario, turn
    
    # Identificamos dinámicamente al otro proceso (0 -> 1, 1 -> 0)
    id_rival = 1 - mi_id  
    
    # 1. PROTOCOLO DE ENTRADA (Algoritmo Dekker)
    flag[mi_id] = True
    
    while flag[id_rival]:
        # Si ambos levantamos la bandera, el turno decide quién cede el paso
        if turn == id_rival:
            flag[mi_id] = False  # Bajo mi bandera por cortesía para evitar Deadlock
            
            while turn == id_rival:
                pass  # Espera activa: me quedo en pausa hasta que sea mi turno
            
            flag[mi_id] = True   # Vuelvo a intentar entrar
            
    # 2. SECCIÓN CRÍTICA (Acceso exclusivo)
    valor_actual = inventario
    
    # Pausa artificial para forzar un error si el algoritmo fallara
    time.sleep(0.1) 
    
    if mi_id == 0:
        inventario = valor_actual + 10
    else:
        inventario = valor_actual - 5
        
    # 3. PROTOCOLO DE SALIDA
    turn = id_rival      # Le paso el turno de desempate al otro proceso
    flag[mi_id] = False  # Anuncio que ya salí de la zona crítica


# Ejecución de la simulación
hilo_0 = threading.Thread(target=modificar_inventario, args=(0,))
hilo_1 = threading.Thread(target=modificar_inventario, args=(1,))

# Iniciamos los procesos simultáneamente
hilo_0.start()
hilo_1.start()

# join() obliga al programa principal a esperar que ambos hilos terminen
hilo_0.join()
hilo_1.join()

print(f"Valor final correcto del inventario: {inventario}")