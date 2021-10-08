
using UnityEngine;
using UnityEngine.UI;

public class Inventory : MonoBehaviour
{
    public bool[] isFull;
    public GameObject[] slots;
    GameManager manager;
    public GameObject[] itemTypes;
    
    public void Start()
    {
        manager = GameObject.FindGameObjectWithTag("Game Manager").GetComponent<GameManager>();
    }
    
    public void OnButtonTapped(int slot)
    {
        if (slots[slot].transform.childCount > 0)
        {
            GameObject child = slots[slot].transform.GetChild(0).gameObject;
            if (child != null)
            {
                switch (child.GetComponent<PickUp>().itemType)
                {
                    case PickUp.ItemType.HealthPotion:
                        manager.IncreasePlayerhealth(50);
                        break;
                    case PickUp.ItemType.StaminaPotion:
                        break;
                    default:
                        break;
                }

                Destroy(child);
                isFull[slot] = false;
            }
        }
    }
}
